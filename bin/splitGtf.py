# Adapted from https://github.com/ablab/IsoQuant
from enum import Enum, unique

import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="Split GTF file into full, known and novel.")
parser.add_argument('-g','--gtf', required=True, type=str, help="annotation file in GTF format")
parser.add_argument('-t','--tool', required=True, type=str, help="tool used to generate the GTF file")
parser.add_argument('-o','--outdir', type=str, default='output', help="output directory")

class TranscriptType(Enum):
    known = 1
    novel = 2
    undefined = 0

class StringTieSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        #if l.find("reference_id") != -1:
        if l.find("cmp_ref") != -1:
            return TranscriptType.known
        else:
            return TranscriptType.novel

class BambuSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.find('novelTranscript "FALSE"') != -1:
            return TranscriptType.known
        else:
            return TranscriptType.novel
        
class IsoscelesSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.find('compatible_tx') != -1:
            return TranscriptType.known
        else:
            return TranscriptType.novel


class IsoQuantSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.find(".nic") != -1 or l.find(".nnic") != -1:
            return TranscriptType.novel
        elif l.find('transcript_id "SIRV') != -1 or l.find('transcript_id "ENS') != -1:
            return TranscriptType.known
        return TranscriptType.undefined

class FlamesSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.split('\t')[1] == "reference":
            return TranscriptType.known
        else:
            return TranscriptType.novel


SEPARATE_FUNCTORS = {'sockeye':StringTieSeparator,
                     'stringtie':StringTieSeparator,
                    'flames':FlamesSeparator,
                    'bambu':BambuSeparator,
                    'isoquant':IsoQuantSeparator,
                    'isoquant-sockeye':IsoscelesSeparator,
                    'isosceles':IsoscelesSeparator,
                    'isosceles-sockeye':IsoscelesSeparator}

def split_gtf(ingtf_path, seaprator, out_full_path, out_known_path, out_novel_path):
    out_full = open(out_full_path, "w")
    out_known = open(out_known_path, "w")
    out_novel = open(out_novel_path, "w")

    open_fn = gzip.open if ingtf_path.endswith(".gz") else open
    for l in open_fn(ingtf_path, "rt"):
        if l.startswith("#"):
            continue
        ttype = seaprator.separate(l)
        if ttype == TranscriptType.undefined:
            continue
        out_full.write(l)
        if ttype == TranscriptType.novel:
            out_novel.write(l)
        elif ttype == TranscriptType.known:
            out_known.write(l)
    out_full.close()
    out_novel.close()
    out_known.close()

def split_gtf2(ingtf_path, seaprator, out_full_path, out_known_path, out_novel_path):
    out_full = open(out_full_path, "w")
    out_known = open(out_known_path, "w")
    out_novel = open(out_novel_path, "w")

    open_fn = gzip.open if ingtf_path.endswith(".gz") else open
    transcript_lines = {}
    for l in open_fn(ingtf_path, "rt"):
        if l.startswith("#"):
            continue
        fields = l.split('\t')
        if fields[2] == "transcript":
            ttype = seaprator.separate(l)
            if ttype == TranscriptType.undefined:
                continue
            transcript_id = [x for x in fields[8].split(';') if 'transcript_id' in x][0].split('"')[1]
            transcript_lines[transcript_id] = (ttype, [l])
        elif fields[2] == "exon":
            transcript_id = [x for x in fields[8].split(';') if 'transcript_id' in x][0].split('"')[1]
            if transcript_id in transcript_lines:
                transcript_lines[transcript_id][1].append(l)

    for transcript_id, (ttype, lines) in transcript_lines.items():
        for line in lines:
            out_full.write(line)
            if ttype == TranscriptType.novel:
                out_novel.write(line)
            elif ttype == TranscriptType.known:
                out_known.write(line)

    out_full.close()
    out_novel.close()
    out_known.close()

def main():
    args = parser.parse_args()
    gtf_path = args.gtf
    tool = args.tool
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    out_full_path = os.path.join(outdir, "full.gtf")
    out_known_path = os.path.join(outdir, "known.gtf")
    out_novel_path = os.path.join(outdir, "novel.gtf")
    seaprator = SEPARATE_FUNCTORS[tool](gtf_path)

    if tool == 'sockeye':
        split_gtf2(gtf_path, seaprator, out_full_path, out_known_path, out_novel_path)
    else:
        split_gtf(gtf_path, seaprator, out_full_path, out_known_path, out_novel_path)

if __name__ == "__main__":
    main()