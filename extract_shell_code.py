#!/usr/bin/env python3

import re
from pathlib import Path
import argparse
import logging as log

def extract_shell_code(markdown_file: Path, output_file: Path):
    with markdown_file.open('r', encoding='utf-8') as md_file:
        content = md_file.read()
        
    # Regex to match ```shell ... ``` and ``` shell ... ```
    shell_code_blocks = re.findall(r'```[ ]?shell\n(.*?)```', content, re.DOTALL)
    
    with output_file.open('w', encoding='utf-8') as bash_file:
        bash_file.write('#!/bin/bash\n\n')
        for block in shell_code_blocks:
            bash_file.write(block.strip() + '\n')
    
    # Make the output file executable
    output_file.chmod(output_file.stat().st_mode | 0o111)

def main():
    parser = argparse.ArgumentParser(description="Extract shell code blocks from a Markdown file and save to an executable bash file.")
    parser.add_argument('markdown_file', type=Path, help="Path to the input Markdown file.")
    parser.add_argument('output_file', type=Path, help="Path to the output bash file.")
    args = parser.parse_args()
    
    log.basicConfig(level=log.INFO, format='%(levelname)s: %(message)s')
    log.info(f"Extracting shell code from {args.markdown_file} to {args.output_file}")
    
    extract_shell_code(args.markdown_file, args.output_file)
    log.info(f"Shell code extraction complete. Output saved to {args.output_file}")

if __name__ == "__main__":
    main()
