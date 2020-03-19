VERSION 1.2.1
-------------------------
 - Move all scripts into `code` folder
 - Solve bug in `make_graphs.R`: (no limit on max cells)

VERSION 1.2.0
-------------------------
 - Added two additional tags that report base quality of Cell Barcode (XQ) and UMI (Xq)

VERSION 1.1.0
-------------------------
- Removed separate `barcode.txt` file which is now hard-coded into the script

Other stuff
-------------------------
- `split_bam.py`: a script to create separate bam files for every single cell.
  (requires [simplesam](https://github.com/mdshw5/simplesam))
