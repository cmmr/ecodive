// Runs Qurro 0.9.0's own JavaScript (support_files/js/feature_computation.js,
// BSD-3, Qurro development team) under Node.js, for ecodive's
// log_ratio_features() and log_ratio() fixtures. Called by
// data-raw/expected/log_ratio_plots.py, not by hand:
//
//   node log_ratio_plots.js <qurro package dir> <cases.json> <out.json>
//
// For each case of cases.json ({id, n, ranks: [{"Feature ID", r}], counts:
// {feature: {sample: count}}, samples}), writes:
//
//   top, bot   Qurro's autoselection: filterFeatures(rankPlotJSON, n, "r",
//              "autoLiteralTop" / "autoLiteralBot"), as feature IDs in the
//              order Qurro returns them.
//   balance    Each sample's log ratio as Qurro's sample plot computes it:
//              updateBalanceMulti()'s sums over the selected features, in
//              their returned order, then computeBalance(top, bot) (null
//              where a sum is 0).
//
// feature_computation.js is an AMD module that needs only dom_utils'
// getNumberIfValid() (which calls vega.toNumber()); a minimal define() and a
// dom_utils stand-in with the same behaviour for plain numeric strings load
// it without a browser or Vega.

"use strict";
const fs = require("fs");
const path = require("path");

const [qurroDir, casesFile, outFile] = process.argv.slice(2);

const domUtils = {
    getNumberIfValid: function (val) {
        if (typeof val === "string") {
            const t = val.trim();
            const nval = t === "" ? null : +t;
            if (nval !== null && isFinite(nval)) return nval;
        } else if (typeof val === "number") {
            if (isFinite(val)) return val;
        }
        return NaN;
    },
};

let fc = null;
global.define = function (deps, factory) {
    fc = factory(domUtils);
};
const src = fs.readFileSync(
    path.join(qurroDir, "support_files", "js", "feature_computation.js"),
    "utf8"
);
// eslint-disable-next-line no-eval
(0, eval)(src);

const cases = JSON.parse(fs.readFileSync(casesFile, "utf8"));
const out = [];

for (const cs of cases) {
    // The rank plot's dataset, as Qurro's JSON holds it.
    const rankPlotJSON = {
        data: { name: "data" },
        datasets: {
            data: cs.ranks.map((row) => ({ "Feature ID": row.id, r: row.r })),
            qurro_rank_ordering: ["r"],
            qurro_feature_metadata_ordering: [],
        },
    };
    const top = fc.filterFeatures(rankPlotJSON, String(cs.n), "r", "autoLiteralTop");
    const bot = fc.filterFeatures(rankPlotJSON, String(cs.n), "r", "autoLiteralBot");

    // display.js: sumAbundancesForSampleFeatures(), then computeBalance().
    const sum = (feats, s) => {
        let a = 0;
        for (let t = 0; t < feats.length; t++) {
            const f = feats[t]["Feature ID"];
            const c = cs.counts[f][s];
            a += c === undefined ? 0 : c;
        }
        return a;
    };
    const balance = cs.samples.map((s) => fc.computeBalance(sum(top, s), sum(bot, s)));

    out.push({
        id: cs.id,
        top: top.map((f) => f["Feature ID"]),
        bot: bot.map((f) => f["Feature ID"]),
        balance: balance,
    });
}

fs.writeFileSync(outFile, JSON.stringify(out));
