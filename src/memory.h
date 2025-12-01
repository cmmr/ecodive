// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED

void  init_n_ptrs(int n);
void* safe_malloc(size_t bytes);
int   is_safe_ptr(void *ptr);
void  free_all(void);
void* free_one(void *ptr);
void* maybe_free_one(void *ptr);

#endif
