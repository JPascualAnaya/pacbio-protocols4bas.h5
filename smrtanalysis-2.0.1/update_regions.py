#!/bin/env python
"""
Given a pair of new regions for each hole number that matched vector sequence,
update the hole number's original region entry with the maximum of the two new
regions.
"""
usage = """%prog filename_prefix new_regions.tab"""

from collections import defaultdict
import csv
import h5py
from optparse import OptionParser
import os
import sys


def update_regions(filename_prefix, regions_filename):
    """
    Given a tab-delimited file of regions per hole, find each hole's original
    high quality region and resize that region with the largest of the two new
    regions when applicable.
    """
    filename = filename_prefix

    # Read in region coordinates and group them by filename and then by hole number.
    regions_by_file = defaultdict(dict)

    with open(regions_filename, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            # Each row looks like this with the first column containing the
            # filename/hole_number and the other columns the start and end of
            # the non-vector region.
            # m120301_070511_42141_c100304962550000001523012408061245_s1_p0/10 418 934
            pieces = row[0].split("/")
            smrtcell_name, hole_number = pieces[:2]
            hole_number = int(hole_number)
            row_ints = map(int, row[1:])
            regions_by_file[filename].setdefault(hole_number, []).append(row_ints)

    # Open each unique base file and search regions by hole number.
    for filename, regions_to_update in regions_by_file.iteritems():
        sys.stderr.write("Found %i regions to update for %s.masked.bax.h5\n" % (len(regions_to_update), filename))

        # Use the first and last hole number to reduce the search space of the
        # base file.
        first_hole = min(regions_to_update.keys())
        last_hole = max(regions_to_update.keys())

        # Open the base file in append mode to allow in-place updates to the
        # Regions dataset.
        total_bases_removed = 0
        hdf5_filename = "%s.masked.bax.h5" % filename

        if not os.path.exists(hdf5_filename):
            raise Exception("File '%s' doesn't exist." % hdf5_filename)

        with h5py.File(hdf5_filename, "a") as h5f:
            for i in xrange(len(h5f["PulseData/Regions"])):
                region = h5f["PulseData/Regions"][i]
                hole_number = region[0]

                if hole_number < first_hole or region[2] == region[3]:
                    # Skip regions that are either before the first hole we want
                    # or that have no high quality region.
                    continue
                elif hole_number > last_hole:
                    # Stop searching Regions once we've passed the last hole in
                    # our search set.
                    break
                elif hole_number in regions_to_update:
                    # Find all new regions that intersect this hole's high
                    # quality region, keep only those regions where start <= end.
                    new_regions = filter(
                        lambda x: x[0] <= x[1],
                        [intersect_regions(region[2:4], region_to_update)
                         for region_to_update in regions_to_update[hole_number]]
                    )

                    # If there are any regions that passed the above tests, pick
                    # the region with the largest range (based on end -
                    # start). Otherwise, there are no valid new regions within
                    # this high quality region, so the entire region must be set
                    # to a zero length.
                    if len(new_regions) > 0:
                        new_region = max(new_regions, key=lambda r: r[1] - r[0])
                    else:
                        new_region = (0, 0)

                    # Don't try to update regions with the same range.
                    if tuple(region[2:4]) != tuple(new_region):
                        sys.stderr.write("update %i: %s to %s\n" % (hole_number, tuple(region[2:4]), tuple(new_region)))
                        total_bases_removed += (region[3] - region[2]) - (new_region[1] - new_region[0])
                        region[2:4] = new_region
                        h5f["PulseData/Regions"][i] = region

            sys.stderr.write("total bases removed: %i\n" % total_bases_removed)


def intersect_regions(region_a, region_b):
    """
    Given a pair of coordinates, returns the coordinates that overlap between
    the two regions.
    """
    return (max(region_a[0], region_b[0]), min(region_a[1], region_b[1]))


if __name__ == "__main__":
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if len(args) != 2:
        sys.stderr.write("Specify a filename prefix and a tab-delimited regions file.\n")
        parser.print_usage()
        sys.exit(1)

    update_regions(*args)
