// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

#include "ecodive.h"
#include <limits.h>   // INT_MAX
#include <stdarg.h>   // va_list
#include <stdio.h>    // vsnprintf

#ifndef NORET
#  define NORET
#endif


/*
 * Newick parser
 * =============
 *
 * The string is scanned twice, front to back, so run time is linear in its
 * length regardless of the tree's shape:
 *
 *   1. count_tree() counts '(' (internal nodes) and ',' (leaves - 1), and
 *      checks that quotes, comments, and parentheses are balanced.
 *
 *   2. C_read_tree() walks the string with an explicit stack (no recursion,
 *      so arbitrarily deep trees are fine), filling in the edge matrix,
 *      edge lengths, and labels.
 *
 * Nodes, leaves, and edges are numbered in the order they are encountered,
 * i.e. pre-order / "cladewise", the same as ape::read.tree().
 *
 * Lexical rules:
 *
 *   - A label may be enclosed in single quotes. Inside quotes every byte is
 *     literal - including ( ) [ ] : ; , ' _ blanks and newlines - except that
 *     two adjacent single quotes ('') represent one literal apostrophe.
 *
 *   - Outside of quotes:
 *       * [...] is a comment and is ignored. Comments may nest, and quotes
 *         inside a comment have no special meaning.
 *       * Blanks around labels and branch lengths are ignored.
 *       * A line break, together with any blanks touching it, is removed
 *         entirely, so a label or number wrapped onto the next line is
 *         rejoined (e.g. "0.1\n  2" is read as "0.12").
 *       * Underscores in labels become spaces, unless `underscores` is TRUE.
 *
 *   - Only the first tree (up to the first ';') is read. A length on the
 *     root edge is accepted but discarded.
 */


typedef struct {
  const char *str;  // the Newick string (not NUL-terminated in our usage)
  size_t      n;    // its length in bytes
  size_t      pos;  // current read position
  cetype_t    ce;   // encoding of str; also used for the output labels
  char       *buf;  // scratch space for the current label or branch length
  size_t      cap;  // capacity of buf, not counting a terminating NUL
  size_t      len;  // number of bytes currently in buf
} parser_t;


#define NO_EDGE ((size_t) -1)

typedef struct {
  size_t node;  // 0-based index of an open internal node
  size_t edge;  // row of the edge leading into it (NO_EDGE for the root)
} frame_t;


NORET static void parse_error(const parser_t *p, size_t pos, const char *fmt, ...) {

  char msg[256];
  va_list args;
  va_start(args, fmt);
  vsnprintf(msg, sizeof(msg), fmt, args);
  va_end(args);

  // Report a 1-based character position rather than a byte offset.
  int    utf8 = p->ce != CE_LATIN1 && p->ce != CE_BYTES;
  double chr  = 1;
  for (size_t i = 0; i < pos && i < p->n; i++)
    if (!utf8 || ((unsigned char) p->str[i] & 0xC0) != 0x80) chr++;

  Rf_error("Invalid Newick string: %s (at character %.0f).", msg, chr);
}


static inline int is_blank   (char c) { return c == ' '  || c == '\t'; }
static inline int is_newline (char c) { return c == '\n' || c == '\r'; }


// Append one byte to the label buffer. count_tree() sizes the buffer to
// the longest possible token, so this check should never fail.
static inline void buf_push(parser_t *p, char c) {
  if (p->len >= p->cap)
    parse_error(p, p->pos, "internal parser error (label buffer overflow)"); // # nocov
  p->buf[p->len++] = c;
}


// p->pos is at a '['. Advance past the matching ']'.
static void skip_comment(parser_t *p) {

  size_t start = p->pos;
  size_t depth = 0;

  do {
    if (p->pos >= p->n)
      parse_error(p, start, "'[' comment is never closed with ']'");
    char c = p->str[p->pos++];
    if      (c == '[') depth++;
    else if (c == ']') depth--;
  } while (depth > 0);
}


// Advance past any blanks, newlines, and comments.
static void skip_ignorable(parser_t *p) {

  while (p->pos < p->n) {
    char c = p->str[p->pos];
    if      (is_blank(c) || is_newline(c)) p->pos++;
    else if (c == '[')                     skip_comment(p);
    else                                   break;
  }
}


// p->pos is at an opening single quote. Advance past the closing quote,
// decoding the label into p->buf when `store` is true.
static void read_quoted(parser_t *p, int store) {

  size_t start = p->pos++;
  p->len = 0;

  for (;;) {

    if (p->pos >= p->n)
      parse_error(
        p, start,
        "quoted label is never closed (write a literal apostrophe "
        "as '' inside a quoted label)" );

    char c = p->str[p->pos++];

    if (c == '\'') {
      if (p->pos < p->n && p->str[p->pos] == '\'') p->pos++; // '' ==> '
      else                                         return;   // closing quote
    }

    if (store) buf_push(p, c);
  }
}


// Read an unquoted label or branch length into p->buf, stopping at the next
// Newick punctuation mark. Comments are dropped, line breaks (and blanks
// touching them) are removed, and trailing blanks are trimmed. The caller
// has already skipped any leading blanks.
static void read_unquoted(parser_t *p) {

  p->len = 0;

  while (p->pos < p->n) {

    char c = p->str[p->pos];

    if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '\'')
      break;

    if (c == '[') {
      skip_comment(p);
    }
    else if (c == ']') {
      // Already rejected by count_tree(), which scans the same text first.
      parse_error(p, p->pos, "']' without a matching '['"); // # nocov
    }
    else if (is_newline(c)) {
      while (p->len > 0 && is_blank(p->buf[p->len - 1])) p->len--;
      while (p->pos < p->n && (is_blank(p->str[p->pos]) || is_newline(p->str[p->pos])))
        p->pos++;
    }
    else {
      buf_push(p, c);
      p->pos++;
    }
  }

  while (p->len > 0 && is_blank(p->buf[p->len - 1])) p->len--;
}


// Read an optional label into p->buf (p->len is 0 when there's no label).
static void read_label(parser_t *p, int underscores) {

  skip_ignorable(p);

  if (p->pos < p->n && p->str[p->pos] == '\'') {
    read_quoted(p, 1); // Quoted labels are never altered.
  }
  else {
    read_unquoted(p);
    if (!underscores)
      for (size_t i = 0; i < p->len; i++)
        if (p->buf[i] == '_') p->buf[i] = ' ';
  }
}


// Read an optional ":length". Returns 0 when there is no branch length.
static double read_length(parser_t *p) {

  skip_ignorable(p);
  if (p->pos >= p->n || p->str[p->pos] != ':') return 0;

  size_t start = p->pos++;
  skip_ignorable(p);
  read_unquoted(p);

  if (p->len == 0)
    parse_error(p, start, "missing branch length after ':'");

  p->buf[p->len] = '\0';
  char  *end;
  double x = strtod(p->buf, &end);

  if (end != p->buf + p->len)
    parse_error(p, start + 1, "invalid branch length '%.40s'", p->buf);

  return x;
}


// Pass 1: count internal nodes and leaves in the first tree, and make sure
// that quotes, comments, and parentheses are all properly closed. Also
// finds the deepest nesting of parentheses, and the longest stretch of text
// between two of ( ) , : ; since no label or branch length can be longer.
static void count_tree(
    parser_t *p, size_t *n_nodes, size_t *n_leafs,
    size_t *max_depth, size_t *max_token ) {

  size_t depth = 0, n_open = 0, n_comma = 0, deepest = 0;
  size_t seg_start = p->pos, max_seg = 0;

  while (p->pos < p->n) {

    char c = p->str[p->pos];

    if (c == '\'') { read_quoted(p, 0); continue; }
    if (c == '[')  { skip_comment(p);   continue; }
    if (c == ']')  parse_error(p, p->pos, "']' without a matching '['");
    if (c == ';')  break;

    if (c == '(' || c == ')' || c == ',' || c == ':') {
      if (p->pos - seg_start > max_seg) max_seg = p->pos - seg_start;
      seg_start = p->pos + 1;
    }

    if (c == '(') {
      n_open++;
      if (++depth > deepest) deepest = depth;
    }
    if (c == ',') n_comma++;
    if (c == ')') {
      if (depth == 0)
        parse_error(p, p->pos, "unbalanced parentheses; ')' without a matching '('");
      depth--;
    }

    p->pos++;
  }

  if (depth > 0)
    parse_error(
      p, p->pos, "unbalanced parentheses; %.0f '(' never closed with ')'",
      (double) depth );

  if (p->pos - seg_start > max_seg) max_seg = p->pos - seg_start;

  *n_nodes   = n_open;
  *n_leafs   = n_comma + 1;
  *max_depth = deepest;
  *max_token = max_seg;
}


SEXP C_read_tree(SEXP sexp_tree_str, SEXP sexp_underscores) {

  if (!isString(sexp_tree_str) || XLENGTH(sexp_tree_str) != 1 ||
      STRING_ELT(sexp_tree_str, 0) == NA_STRING)
    Rf_error("Newick tree must be a single, non-NA character string."); // # nocov

  SEXP sexp_tree_chr = STRING_ELT(sexp_tree_str, 0);
  int  underscores   = asLogical(sexp_underscores) == TRUE;

  parser_t p;
  p.str = CHAR(sexp_tree_chr);
  p.n   = (size_t) LENGTH(sexp_tree_chr);
  p.pos = 0;
  p.ce  = getCharCE(sexp_tree_chr);
  p.buf = NULL;
  p.cap = 0;
  p.len = 0;


  // The tree proper must begin with '(', optionally preceded by comments.
  skip_ignorable(&p);
  if (p.pos >= p.n)
    parse_error(&p, p.pos, "no tree found");
  if (p.str[p.pos] != '(')
    parse_error(&p, p.pos, "a tree must begin with '('");


  // Pass 1: size the output.
  size_t root_pos = p.pos;
  size_t n_nodes, n_leafs, max_depth, max_token;
  count_tree(&p, &n_nodes, &n_leafs, &max_depth, &max_token);

  // R strings are at most INT_MAX bytes, which is too short to trigger this.
  if (n_nodes + n_leafs > INT_MAX)
    parse_error(&p, root_pos, "too many nodes"); // # nocov

  size_t n_edges = n_nodes + n_leafs - 1;


  // Prepare a phylo-compatible data structure to return to R.
  SEXP sexp_result          = PROTECT(allocVector(VECSXP, 5));
  SEXP sexp_edge_mtx        = PROTECT(allocMatrix(INTSXP,  (int) n_edges, 2));
  SEXP sexp_edge_length_vec = PROTECT(allocVector(REALSXP, (R_xlen_t) n_edges));
  SEXP sexp_node_label_vec  = PROTECT(allocVector(STRSXP,  (R_xlen_t) n_nodes));
  SEXP sexp_leaf_label_vec  = PROTECT(allocVector(STRSXP,  (R_xlen_t) n_leafs));

  SET_VECTOR_ELT(sexp_result, 0, sexp_edge_mtx);
  SET_VECTOR_ELT(sexp_result, 1, ScalarInteger((int) n_nodes));
  SET_VECTOR_ELT(sexp_result, 2, sexp_leaf_label_vec);
  SET_VECTOR_ELT(sexp_result, 3, sexp_edge_length_vec);
  SET_VECTOR_ELT(sexp_result, 4, sexp_node_label_vec);

  // Labels default to "". Only non-empty labels are written below.
  for (size_t i = 0; i < n_nodes; i++) SET_STRING_ELT(sexp_node_label_vec, i, R_BlankString);
  for (size_t i = 0; i < n_leafs; i++) SET_STRING_ELT(sexp_leaf_label_vec, i, R_BlankString);

  int    *edge_mtx        = INTEGER(sexp_edge_mtx);
  double *edge_length_vec = REAL(sexp_edge_length_vec);
  memset(edge_mtx,        0, n_edges * 2 * sizeof(int));
  memset(edge_length_vec, 0, n_edges * sizeof(double));

  // Scratch memory; R_alloc is released automatically, even on error.
  p.cap          = max_token;
  p.buf          = R_alloc(max_token + 1, sizeof(char));
  frame_t *stack = (frame_t *) R_alloc(max_depth, sizeof(frame_t));


  // Pass 2: build the tree.
  size_t next_node = 0, next_leaf = 0, next_edge = 0, depth = 0;
  int    expect_child = 1;

  p.pos = root_pos;

  for (;;) {

    skip_ignorable(&p);
    char c = (p.pos < p.n) ? p.str[p.pos] : '\0';


    if (expect_child) {

      if (c == '(') {
        // Open a new internal node.

        if (next_node >= n_nodes || depth >= max_depth ||
            (depth > 0 && next_edge >= n_edges))
          parse_error(&p, p.pos, "internal parser error"); // # nocov

        size_t node = next_node++;
        size_t edge = NO_EDGE;

        if (depth > 0) {
          edge = next_edge++;
          edge_mtx[edge]           = (int) (n_leafs + stack[depth - 1].node + 1);
          edge_mtx[edge + n_edges] = (int) (n_leafs + node + 1);
        }

        stack[depth].node = node;
        stack[depth].edge = edge;
        depth++;
        p.pos++;

        continue; // The first child of this node comes next.
      }

      // Anything else is a leaf, possibly with an empty label.

      if (next_leaf >= n_leafs || next_edge >= n_edges)
        parse_error(&p, p.pos, "internal parser error"); // # nocov

      size_t leaf = next_leaf++;
      size_t edge = next_edge++;
      edge_mtx[edge]           = (int) (n_leafs + stack[depth - 1].node + 1);
      edge_mtx[edge + n_edges] = (int) (leaf + 1);

      read_label(&p, underscores);
      if (p.len > 0)
        SET_STRING_ELT(sexp_leaf_label_vec, leaf, mkCharLenCE(p.buf, (int) p.len, p.ce));

      edge_length_vec[edge] = read_length(&p);

      expect_child = 0;
      continue;
    }


    if (c == ',') {
      // Another child of the current internal node follows.
      p.pos++;
      expect_child = 1;
      continue;
    }


    if (c == ')') {
      // Close the current internal node, then read its label and length.

      if (depth == 0)
        parse_error(&p, p.pos, "internal parser error"); // # nocov

      depth--;
      size_t node = stack[depth].node;
      size_t edge = stack[depth].edge;
      p.pos++;

      read_label(&p, underscores);
      if (p.len > 0)
        SET_STRING_ELT(sexp_node_label_vec, node, mkCharLenCE(p.buf, (int) p.len, p.ce));

      double length = read_length(&p);
      if (edge != NO_EDGE) edge_length_vec[edge] = length;

      if (depth == 0) break; // Closed the root.
      continue;
    }


    // Unreachable for input that passed count_tree(), but be defensive.
    if (p.pos >= p.n || c == ';')
      parse_error(&p, p.pos, "unexpected end of tree"); // # nocov

    if (c == '\'')
      parse_error(&p, p.pos, "unexpected quote; quotes must enclose an entire label");
    if (c == ':')
      parse_error(&p, p.pos, "unexpected ':'; only one branch length is allowed per node");
    if (c == '(')
      parse_error(&p, p.pos, "unexpected '('; expected ',' or ')' first");

    parse_error(&p, p.pos, "unexpected text after a label; expected ',' or ')'");
  }


  // Nothing but comments/whitespace may separate the root from ';'.
  skip_ignorable(&p);
  if (p.pos < p.n && p.str[p.pos] != ';') {
    if (p.str[p.pos] == ',')
      parse_error(
        &p, p.pos,
        "unexpected ',' after the root node; is the tree missing its outermost parentheses?" );
    parse_error(&p, p.pos, "unexpected text after the root node; expected ';'");
  }

  UNPROTECT(5);
  return sexp_result;
}
