#!/usr/bin/env python3
"""
Stacked 1D histograms: data on the stacked backgrounds, signals as lines, and a data/MC ratio.

    python3 python/stack.py METanalysis_TestPUppiMETCut_Eta2p4_Fpix
    python3 python/stack.py 'METanalysis_TestPUppiMETCut_Eta2p4_*' --logy --rebin 2 --signal-scale 100
    python3 python/stack.py --list 'Nm1_event_*'                     # the 1D histograms matching, nothing drawn
    python3 python/stack.py EventCutflow --bkg qcd,ttbar --data '' --signal gluino   # no data

Histograms are given by name or fnmatch pattern (quote patterns), matched against the 1D histograms of all
the samples (those of the data file first). The samples are inputs of
cfg/plots.json (--config): its "stack" block, or --data NAME, --bkg A,B,... (stacking order, the first at
the bottom) and --signal A,B,... on the command line; labels, colours and the luminosity text come from the
same file. A sample without the histogram is left out of that plot, with a warning.

One file per histogram in <output_dir>/stack/ (--output-dir), named after the histogram, with _log for
--logy. The integral of each sample is printed. The ratio panel shows data / sum of the backgrounds, with
the statistical uncertainty of the backgrounds as a band; it is left out with --no-ratio, or without data
or backgrounds. Needs PyROOT: source env.sh first.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hscp import plots, stack  # noqa: E402


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hists", nargs="*", help="histogram names or fnmatch patterns")
    ap.add_argument("--config", default=plots.PLOTS_FILE, help="inputs, labels, colours (default: %(default)s)")
    ap.add_argument("--data", help="input drawn as data ('' for none; default: the stack block)")
    ap.add_argument("--bkg", help="comma-separated inputs stacked, the first at the bottom")
    ap.add_argument("--signal", help="comma-separated inputs drawn as lines")
    ap.add_argument("--signal-scale", type=float, default=1.0, help="factor applied to the signals (shown in the legend)")
    ap.add_argument("--logy", action="store_true", help="logarithmic y axis")
    ap.add_argument("--rebin", type=int, default=1, help="merge N bins into one")
    ap.add_argument("--xmin", type=float)
    ap.add_argument("--xmax", type=float)
    ap.add_argument("--ymin", type=float, help="y minimum (default: 0, or from the contents with --logy)")
    ap.add_argument("--ymax", type=float, help="y maximum (default: from the contents)")
    ap.add_argument("--xtitle", help="x-axis title (default: the one of the histogram, else its name)")
    ap.add_argument("--no-ratio", action="store_true", help="no data/MC panel")
    ap.add_argument("--output-dir", help="default: <output_dir of the config>/stack")
    ap.add_argument("--format", default="pdf", help="comma-separated file types, e.g. pdf,png (default: pdf)")
    ap.add_argument("--list", action="store_true", help="list the 1D histograms matching, draw nothing")
    return ap.parse_args(argv)


class Drawer:
    def __init__(self, ROOT, cfg, samples, files, args):
        self.R, self.cfg, self.s, self.files, self.a = ROOT, cfg, samples, files, args
        self.bkg_colors = self._root_colors(stack.colors(cfg, samples.backgrounds, stack.BACKGROUND_COLORS))
        self.sig_colors = self._root_colors(stack.colors(cfg, samples.signals, stack.SIGNAL_COLORS))
        self.keep = []                                   # PyROOT objects drawn on the canvas, kept alive

    def _root_colors(self, specs):
        out = {}
        for n, c in specs.items():
            if isinstance(c, int):
                out[n] = c
                continue
            base = getattr(self.R, c[0], None)
            if not isinstance(base, int):
                raise plots.PlotError(f"input {n}: unknown ROOT colour {c[0]}")
            out[n] = base + c[1]
        return out

    def get(self, sample, name):
        h = self.files[sample].Get(name)
        if not h or not h.InheritsFrom("TH1") or h.GetDimension() != 1:
            return None
        h = h.Clone(f"{name}__{sample}")
        h.SetDirectory(0)
        if self.a.rebin > 1:
            h.Rebin(self.a.rebin)
        self.keep.append(h)
        return h

    def draw(self, name, out_base):
        R, a, s = self.R, self.a, self.s
        data = self.get(s.data, name) if s.data else None
        bkgs = [(n, self.get(n, name)) for n in s.backgrounds]
        sigs = [(n, self.get(n, name)) for n in s.signals]
        missing = [n for n, h in bkgs + sigs if h is None] + ([s.data] if s.data and data is None else [])
        bkgs = [(n, h) for n, h in bkgs if h]
        sigs = [(n, h) for n, h in sigs if h]
        if missing:
            print(f"warning: {name} is not a 1D histogram of {', '.join(missing)}: left out")
        drawn = ([data] if data else []) + [h for _, h in bkgs + sigs]
        if not drawn:
            return None

        # the stack as cumulative sums (background k on top of 0..k-1), drawn tallest first: the same picture as
        # a THStack, which crashes in THStack::BuildAndPaint with several histograms in some ROOT versions (6.40)
        total, cumulative = None, []
        for n, h in bkgs:
            h.SetFillColor(self.bkg_colors[n])
            h.SetLineColor(self.bkg_colors[n])
            total = h.Clone(f"{name}__cum_{n}") if total is None else total.Clone(f"{name}__cum_{n}")
            if len(cumulative):
                total.Add(h)
            total.SetDirectory(0)
            total.SetFillColor(self.bkg_colors[n])
            total.SetLineColor(self.bkg_colors[n])
            cumulative.append(total)
        self.keep += cumulative
        for n, h in sigs:
            if a.signal_scale != 1:
                h.Scale(a.signal_scale)
            h.SetLineColor(self.sig_colors[n])
            h.SetLineWidth(3)
            h.SetFillStyle(0)
        if data:
            data.SetMarkerStyle(20)
            data.SetMarkerSize(0.9)
            data.SetLineColor(R.kBlack)
            data.SetMarkerColor(R.kBlack)

        ratio = not a.no_ratio and data is not None and total is not None
        c = R.TCanvas(f"c_{name}", "", 800, 800 if ratio else 650)
        self.keep.append(c)
        if ratio:
            top = R.TPad(f"top_{name}", "", 0, 0.3, 1, 1)
            bottom = R.TPad(f"bottom_{name}", "", 0, 0, 1, 0.3)
            top.SetBottomMargin(0.02)
            bottom.SetTopMargin(0.03)
            bottom.SetBottomMargin(0.35)
            for p in (top, bottom):
                p.SetLeftMargin(0.13)
                p.SetRightMargin(0.04)
                c.cd()
                p.Draw()
            self.keep += [top, bottom]
        else:
            top, bottom = c, None
            c.SetLeftMargin(0.13)
            c.SetRightMargin(0.04)
            c.SetBottomMargin(0.12)
        top.cd()
        top.SetLogy(a.logy)
        top.SetTopMargin(0.08)

        # axis frame, then the stack, its uncertainty, the signals and the data
        frame = drawn[0].Clone(f"{name}__frame")
        frame.SetDirectory(0)
        frame.Reset()
        frame.SetTitle("")
        self.keep.append(frame)
        tops = [h.GetBinContent(i) + (h.GetBinError(i) if h is data else 0)
                for h in ([data] if data else []) + ([total] if total else []) + [h for _, h in sigs]
                for i in range(1, h.GetNbinsX() + 1)]
        positive = [h.GetBinContent(i) for h in drawn for i in range(1, h.GetNbinsX() + 1) if h.GetBinContent(i) > 0]
        hmax = max(tops) if tops else 1.0
        if a.logy:
            ymin = a.ymin if a.ymin is not None else (0.5 * min(positive) if positive else 0.1)
            ymax = a.ymax if a.ymax is not None else (ymin * (hmax / ymin) ** 1.8 if hmax > ymin else 10 * ymin)
        else:
            ymin = a.ymin if a.ymin is not None else 0.0
            ymax = a.ymax if a.ymax is not None else 1.45 * hmax
        frame.SetMinimum(ymin)
        frame.SetMaximum(ymax)
        if a.xmin is not None or a.xmax is not None:
            ax = frame.GetXaxis()
            frame.GetXaxis().SetRangeUser(ax.GetXmin() if a.xmin is None else a.xmin, ax.GetXmax() if a.xmax is None else a.xmax)
        xtitle = a.xtitle or drawn[0].GetXaxis().GetTitle() or name
        frame.GetXaxis().SetTitle(xtitle)
        frame.GetYaxis().SetTitle("Events")
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.045)
        frame.GetYaxis().SetTitleOffset(1.2)
        if ratio:
            frame.GetXaxis().SetLabelSize(0)
            frame.GetXaxis().SetTitleSize(0)
        else:
            frame.GetXaxis().SetTitleSize(0.045)
            frame.GetXaxis().SetLabelSize(0.04)
        frame.Draw("AXIS")
        if bkgs:
            for h in reversed(cumulative):                  # the full sum first, then the lower layers
                h.Draw("HIST SAME")
            band = total.Clone(f"{name}__band")
            band.SetDirectory(0)
            band.SetFillStyle(3354)
            band.SetFillColor(R.kGray + 2)
            band.SetLineColor(R.kGray + 2)
            band.SetMarkerSize(0)
            band.Draw("E2 SAME")
            self.keep.append(band)
        for _, h in sigs:
            h.Draw("HIST SAME")
        if data:
            data.Draw("E1 X0 SAME")
        frame.Draw("AXIS SAME")

        n_entries = (1 if data else 0) + len(bkgs) + len(sigs) + (1 if bkgs else 0)
        leg = R.TLegend(0.55, max(0.45, 0.89 - 0.055 * n_entries), 0.94, 0.89)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        if data:
            leg.AddEntry(data, stack.legend_label(self.cfg, s.data), "pe")
        for n, h in reversed(bkgs):                     # top of the stack first
            leg.AddEntry(h, stack.legend_label(self.cfg, n), "f")
        if bkgs:
            leg.AddEntry(band, "Stat. unc.", "f")
        for n, h in sigs:
            leg.AddEntry(h, stack.legend_label(self.cfg, n, a.signal_scale), "l")
        leg.Draw()
        self.keep.append(leg)

        cms = R.TLatex()
        cms.SetNDC()
        cms.SetTextFont(42)
        cms.SetTextSize(0.045)
        cms.DrawLatex(top.GetLeftMargin(), 1 - top.GetTopMargin() + 0.015, "#bf{CMS} #it{Private work}")
        if self.cfg.lumi_label:
            cms.SetTextAlign(31)
            cms.DrawLatex(1 - top.GetRightMargin(), 1 - top.GetTopMargin() + 0.015, self.cfg.lumi_label)
        self.keep.append(cms)

        if ratio:
            bottom.cd()
            rframe = frame.Clone(f"{name}__rframe")
            rframe.SetDirectory(0)
            rframe.SetMinimum(0.0)
            rframe.SetMaximum(2.0)
            for ax, size in ((rframe.GetXaxis(), 0.11), (rframe.GetYaxis(), 0.1)):
                ax.SetTitleSize(size)
                ax.SetLabelSize(0.09)
            rframe.GetXaxis().SetTitle(xtitle)
            rframe.GetXaxis().SetTitleOffset(1.2)
            rframe.GetYaxis().SetTitle("Data / MC")
            rframe.GetYaxis().SetTitleOffset(0.55)
            rframe.GetYaxis().SetNdivisions(505)
            rframe.Draw("AXIS")
            rband = total.Clone(f"{name}__rband")
            rband.SetDirectory(0)
            r = data.Clone(f"{name}__ratio")
            r.SetDirectory(0)
            for i in range(0, total.GetNbinsX() + 2):
                mc = total.GetBinContent(i)
                rband.SetBinContent(i, 1.0 if mc > 0 else 0.0)
                rband.SetBinError(i, total.GetBinError(i) / mc if mc > 0 else 0.0)
                r.SetBinContent(i, data.GetBinContent(i) / mc if mc > 0 else 0.0)
                r.SetBinError(i, data.GetBinError(i) / mc if mc > 0 else 0.0)
            rband.SetFillStyle(3354)
            rband.SetFillColor(R.kGray + 2)
            rband.SetLineColor(R.kGray + 2)
            rband.SetMarkerSize(0)
            rband.Draw("E2 SAME")
            ax = rframe.GetXaxis()
            line = R.TLine(ax.GetBinLowEdge(ax.GetFirst()), 1.0, ax.GetBinUpEdge(ax.GetLast()), 1.0)
            line.SetLineStyle(2)
            line.Draw()
            r.Draw("E1 X0 SAME")
            self.keep += [rframe, rband, r, line]

        written = []
        for fmt in [f.strip() for f in a.format.split(",") if f.strip()]:
            path = f"{out_base}.{fmt}"
            c.SaveAs(path)
            written.append(path)

        yields = ([(s.data, data.Integral())] if data else []) + [(n, h.Integral()) for n, h in bkgs] \
            + ([("total MC", total.Integral())] if total else []) \
            + [(n + (f" (x{a.signal_scale:g})" if a.signal_scale != 1 else ""), h.Integral()) for n, h in sigs]
        print(f"{name}:  " + ",  ".join(f"{n} {v:.6g}" for n, v in yields))
        return written


def one_d_keys(ROOT, f):
    """The names of the 1D histograms of a file (TProfiles excluded), in file order."""
    out = []
    for k in f.GetListOfKeys():
        cl = ROOT.TClass.GetClass(k.GetClassName())
        if cl and cl.InheritsFrom("TH1") and not cl.InheritsFrom("TH2") and not cl.InheritsFrom("TH3") \
                and not cl.InheritsFrom("TProfile") and k.GetName() not in out:
            out.append(k.GetName())
    return out


def main(argv=None):
    args = parse_args(argv)
    if not args.hists:
        sys.exit("stack: give histogram names or patterns (--list '*' shows them)")
    if args.rebin < 1:
        sys.exit("stack: --rebin must be at least 1")
    try:
        cfg = plots.load(args.config)
        s = stack.samples(cfg, args.data, stack.split_names(args.bkg), stack.split_names(args.signal))
    except plots.PlotError as e:
        sys.exit(f"stack: {e}")
    names = ([s.data] if s.data else []) + s.backgrounds + s.signals
    absent = [f"{n} ({cfg.inputs[n]})" for n in names if not os.path.isfile(cfg.inputs[n])]
    if absent:
        sys.exit("stack: missing input file(s):\n  " + "\n  ".join(absent))

    try:
        import ROOT
    except ImportError:
        sys.exit("stack: PyROOT not found (source env.sh first)")
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    files = {}
    for n in names:
        f = ROOT.TFile.Open(cfg.inputs[n])
        if not f or f.IsZombie():
            sys.exit(f"stack: cannot open {cfg.inputs[n]}")
        files[n] = f

    keys = list(dict.fromkeys(k for n in names for k in one_d_keys(ROOT, files[n])))
    try:
        hists = stack.expand(args.hists, keys)
    except plots.PlotError as e:
        sys.exit(f"stack: {e} among the 1D histograms of {', '.join(names)}")
    if args.list:
        print("\n".join(hists))
        return 0

    out_dir = os.path.abspath(args.output_dir) if args.output_dir else os.path.join(cfg.output_dir, "stack")
    os.makedirs(out_dir, exist_ok=True)
    try:
        drawer = Drawer(ROOT, cfg, s, files, args)
    except plots.PlotError as e:
        sys.exit(f"stack: {e}")
    skipped = 0
    for name in hists:
        base = os.path.join(out_dir, stack.file_name(name) + ("_log" if args.logy else ""))
        if drawer.draw(name, base) is None:
            skipped += 1
            print(f"warning: {name}: no sample holds it, nothing drawn")
    print(f"\n{len(hists) - skipped} plot(s) in {out_dir}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
