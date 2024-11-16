#!/usr/bin/env python3
# Author: Ramin Yousefpour Shahrivar

from collections import Counter
import configparser
import pandas as pd
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def load_config(config_path='config.ini'):
    config_file = configparser.ConfigParser()
    config_file.read(config_path)
    config = {section: dict(config_file.items(section))
              for section in config_file.sections()}
    config['DEFAULT'] = dict(config_file['DEFAULT'])
    return config['DEFAULT']

config = load_config()
units = config['units'].split(' ')
output_dir = Path(config['repository']) / 'results' / \
    config['sample'] / f"composition{config['quality']}"
output_dir.mkdir(parents=True, exist_ok=True)

def count_ribo(file, bg_file, region=None):
    if region:
        region_bed = Path(config['repository']) / 'results' / config['sample'] / \
            f"coordinate{config['quality']}" / \
            f"{config['sample']}-{region}.bed"
        nucs_tab = output_dir / f"{config['sample']}-{region}.nucs.tab"
    else:
        region_bed = Path(config['repository']) / 'results' / config['sample'] / \
            f"coordinate{config['quality']}" / f"{config['sample']}.bed"
        nucs_tab = output_dir / f"{config['sample']}.nucs.tab"

    subprocess.run(
        f"bedtools getfasta -s -fi {config['fasta']} -bed {region_bed} | grep -v '>' > {nucs_tab}",
        shell=True,
        check=True
    )

    ribo_num = Counter(nucs_tab.read_text().replace('\n', '').upper())
    ribo_sum = sum(ribo_num.values())

    bg_freq = pd.read_csv(bg_file, sep='\t', header=None)[1]
    return ribo_sum, bg_freq, dict(ribo_num)

def normalize_ribo(ribo_sum, bg_freq, ribo_num):
    ribo_normal = {
        f"r{base}": round((ribo_num.get(base, 0) / ribo_sum) / bg_freq[idx] * 100 / sum(
            (ribo_num.get(b, 0) / ribo_sum) / bg_freq[idx] for idx, b in enumerate("ACGT")), 5)
        for idx, base in enumerate("ACGT")
    }
    return ribo_normal

def save_files(file, bg_file, region=None):
    ribo_sum, bg_freq, ribo_num = count_ribo(file, bg_file, region)
    norm_ribo = normalize_ribo(ribo_sum, bg_freq, ribo_num)
    if region:
        counts_file = output_dir / f"{config['sample']}-{region}.counts.txt"
        freq_file = output_dir / f"{config['sample']}-{region}.frequencies.txt"
    else:
        counts_file = output_dir / f"{config['sample']}.counts.txt"
        freq_file = output_dir / f"{config['sample']}.frequencies.txt"
    counts_file.write_text(str(ribo_num))
    freq_file.write_text(str(norm_ribo))

def process_region(region=None):
    if region:
        file = output_dir / f"{config['sample']}-{region}.nucs.tab"
        bg_file = Path(
            f"/home/core/ribose-map/{config['fasta']}").stem + f'-{region}.txt'
    else:
        file = output_dir / f"{config['sample']}.nucs.tab"
        bg_file = Path(
            f"/home/core/ribose-map/{config['fasta']}").stem + '.txt'

    save_files(file, bg_file, region)

def main():
    all_regions = units if units else [None]
    with ThreadPoolExecutor() as executor:
        executor.map(process_region, all_regions)

if __name__ == "__main__":
    main()
    print(f"Composition module for {config['sample']} ran successfuly.")
