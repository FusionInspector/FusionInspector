#!/usr/bin/env python3
"""
Test script to demonstrate fusion file preprocessing.
Shows how both simple list and STAR-Fusion output formats are handled.
"""

import tempfile
import os

# Create test files
def test_fusion_preprocessing():

    # Test 1: Simple list format
    simple_list = """# Comment line
ACACA--STAC2
EML4--ALK
BCR--ABL1
"""

    # Test 2: STAR-Fusion output format
    star_fusion_output = """#FusionName	JunctionReadCount	SpanningFragCount	est_J	est_S	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	JunctionReads	SpanningFrags	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots
EML4--ALK	150	50	150	50	ONLY_REF_SPLICE	EML4	chr2:42522656:+	ALK	chr2:29415640:-	read1,read2,read3	frag1,frag2	YES_LDAS	5.23	GT	1.92	AG	1.95	["INTERGENIC"]
ACACA--STAC2	100	30	100	30	ONLY_REF_SPLICE	ACACA	chr17:35518578:-	STAC2	chr17:37839896:+	read4,read5	frag3,frag4	YES_LDAS	3.45	GT	1.85	AG	1.90	["NEIGHBORS"]
EML4--ALK	80	20	80	20	ONLY_REF_SPLICE	EML4	chr2:42522700:+	ALK	chr2:29415640:-	read6,read7	frag5,frag6	YES_LDAS	2.12	GT	1.88	AG	1.93	["INTERGENIC"]
BCR--ABL1	200	60	200	60	ONLY_REF_SPLICE	BCR	chr22:23632601:+	ABL1	chr9:133729450:-	read8,read9	frag7,frag8	YES_LDAS	8.90	GT	1.95	AG	1.98	["INTERCHROMOSOMAL"]
"""

    with tempfile.TemporaryDirectory() as tmpdir:

        # Test simple list
        simple_file = os.path.join(tmpdir, "simple_list.txt")
        with open(simple_file, "w") as f:
            f.write(simple_list)

        print("Test 1: Simple list format")
        print(f"Input file: {simple_file}")
        with open(simple_file) as f:
            print("Contents:")
            print(f.read())
        print()

        # Test STAR-Fusion output
        star_file = os.path.join(tmpdir, "star_fusion_output.tsv")
        with open(star_file, "w") as f:
            f.write(star_fusion_output)

        print("Test 2: STAR-Fusion output format")
        print(f"Input file: {star_file}")
        print("Expected unique fusions to extract: ACACA--STAC2, BCR--ABL1, EML4--ALK (3 unique)")
        print("Note: EML4--ALK appears twice with different breakpoints, should be deduplicated")
        print()

        # Show what the preprocessing would extract
        print("Fusion names that would be extracted:")
        fusions_seen = set()
        with open(star_file) as f:
            first_line = f.readline().strip()
            headers = first_line.split("\t")
            fusion_name_col = headers.index("#FusionName")

            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > fusion_name_col:
                    fusion_name = parts[fusion_name_col]
                    fusions_seen.add(fusion_name)

        for fusion in sorted(fusions_seen):
            print(f"  {fusion}")

        print(f"\nTotal unique fusions: {len(fusions_seen)}")

if __name__ == "__main__":
    test_fusion_preprocessing()
