#!/bin/bash
cd `dirname $0`
# diff  ../../Patch/patchSys1.m patchSys1.m > patchSys1.txt
# diff  ../../Patch/configPatches1.m configPatches1.m > configPatches1.txt
# diff  ../../Patch/patchEdgeInt1.m patchEdgeInt1.m > patchEdgeInt1.txt
# diff  ../../Patch/patchEdgeInt1test.m patchEdgeInt1test.m > patchEdgeInt1test.txt
diff  ../../Patch/patchEdgeInt2.m patchEdgeInt2.m > patchEdgeInt2.txt
diff  ../../Patch/patchEdgeInt2test.m patchEdgeInt2test.m > patchEdgeInt2test.txt
echo finished script
sleep 5
