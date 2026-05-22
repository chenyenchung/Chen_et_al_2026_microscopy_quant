#!/usr/bin/env python
import imagej
import scyjava as sj
from pathlib import Path
import re
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument('--lif-path', required=True)
parser.add_argument('--out-dir', required=True)
parser.add_argument('--channels', nargs='+', type=int, required=True)
parser.add_argument('--suffix', default='', help='If specified, pixel width will not be included in output filename, and this suffix will be appended instead.')
args = parser.parse_args()

# Configure logging to show INFO level messages.
logging.basicConfig(level=logging.INFO, format='%(message)s')
logging.info(f"Exporting channels {args.channels} from {args.lif_path} to {args.out_dir}")

def strip_lif_prefix(title, lif_path):
    """
    Bio-Formats often names opened series like:
        file.lif - stack_name

    We want only:
        stack_name
    """
    title = str(title)
    lif_path = Path(lif_path)

    prefixes = [
        f"{lif_path.name} - ",
        f"{lif_path.stem} - ",
    ]

    for prefix in prefixes:
        if title.startswith(prefix):
            return title[len(prefix):]

    return title


def safe_name(x):
    x = str(x) if x is not None else "unnamed"
    x = re.sub(r"[^\w.\-]+", "_", x).strip("_")
    return x or "unnamed"


def imageplus_array(imps):
    """
    Convert Python list[ImagePlus] -> Java ImagePlus[].
    Needed by RGBStackMerge.mergeChannels().
    """
    arr = sj.jarray(ImagePlus, len(imps))
    for i, imp in enumerate(imps):
        arr[i] = imp
    return arr


def export_lif_channels(lif_path, out_dir, channels=(0,), virtual=True):
    """
    Export selected channels from every series in a LIF file.

    channels=(0,)      -> export first channel as single-channel TIFF
    channels=(1,2,3)   -> export channels 2,3,4 as one 3-channel TIFF

    Channel indices are 0-based, matching Python list indexing.
    """
    lif_path = Path(lif_path).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(channels, int):
        channels = [channels]
    channels = list(channels)

    # First pass: count series.
    reader = ImageReader()
    try:
        reader.setId(str(lif_path))
        n_series = reader.getSeriesCount()
    finally:
        reader.close()

    logging.info(f"{lif_path.name}: {n_series} series")

    for s in range(n_series):
        opts = ImporterOptions()
        opts.setId(str(lif_path))
        opts.setOpenAllSeries(False)
        opts.clearSeries()
        opts.setSeriesOn(s, True)
        opts.setVirtual(virtual)
        opts.setQuiet(True)

        imps = list(BF.openImagePlus(opts))
        if not imps:
            logging.info(f"Series {s}: no ImagePlus returned; skipping")
            continue

        imp = imps[0]
        split_channels = None
        merged = None

        try:
            raw_title = strip_lif_prefix(imp.getTitle(), lif_path)
            name = safe_name(raw_title)

            cal = imp.getCalibration()
            pixel_width = float(cal.pixelWidth)
            unit = cal.getUnit()

            split_channels = list(ChannelSplitter.split(imp))

            n_ch = len(split_channels)
            bad = [c for c in channels if c < 0 or c >= n_ch]
            if bad:
                logging.info(
                    f"Skipping {name}: requested channel(s) {bad}, "
                    f"but series has {n_ch} channel(s)"
                )
                continue

            selected = [split_channels[c] for c in channels]

            if len(selected) == 1:
                out_imp = selected[0]
            else:
                # Merge in exactly the order requested by `channels`.
                merged = RGBStackMerge.mergeChannels(
                    imageplus_array(selected),
                    False,   # do not keep source images duplicated inside ImageJ
                )
                merged.copyScale(imp)
                out_imp = merged

            if args.suffix:
                out_path = out_dir / f"{name}_{args.suffix}.tif"
            else:
                out_path = out_dir / f"{name}_{pixel_width:.2f}.tif"

            # Avoid silent overwrite if stack names collide.
            if out_path.exists():
                ch_tag = "ch" + "-".join(str(c) for c in channels)
                out_path = out_dir / f"{name}_series{s:03d}_{ch_tag}_{pixel_width:.2f}.tif"

            out_imp.setTitle(out_path.stem)
            IJ.saveAsTiff(out_imp, str(out_path))

            logging.info(
                f"Wrote {out_path.name} "
                f"[channels={channels}, pixel_width={pixel_width:.4g} {unit}]"
            )

        finally:
            if merged is not None:
                try:
                    merged.close()
                except Exception:
                    pass

            if split_channels is not None:
                for ch_imp in split_channels:
                    try:
                        ch_imp.close()
                    except Exception:
                        pass

            try:
                imp.close()
            except Exception:
                pass

            sj.gc()

sj.config.add_options('-Xmx6g')
ij = imagej.init('sc.fiji:fiji', mode='headless') # Import required Java classes via scyjava
BF = sj.jimport('loci.plugins.BF')
ImporterOptions = sj.jimport('loci.plugins.in.ImporterOptions')
ImageReader = sj.jimport('loci.formats.ImageReader')
ChannelSplitter = sj.jimport('ij.plugin.ChannelSplitter')
RGBStackMerge = sj.jimport('ij.plugin.RGBStackMerge')
ImagePlus = sj.jimport('ij.ImagePlus')
IJ = sj.jimport('ij.IJ')

export_lif_channels(
    lif_path=args.lif_path,
    out_dir=args.out_dir,
    channels=args.channels,
)