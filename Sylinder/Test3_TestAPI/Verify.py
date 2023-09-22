import sys
import vtk
import glob
import re
import numpy as np
import scipy as sp
import yaml
import copy

# member variables are dynamically added by parsing data files

file = open('RunConfig.yaml')
config = yaml.load(file, Loader=yaml.FullLoader)
file.close()
print(config['endLinkKappa'])
linkKappa = config['endLinkKappa']


class ConBlock(object):
    end0 = None
    end1 = None
    pass


class Frame:

    def __init__(self, conBlockFile=None):
        self.conBlocks = []
        self.parseConBlockFile(conBlockFile)

    def parseFile(self, dataFile, objType, objList):
        # print("Parsing data from " + dataFile)
        # create vtk reader
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(dataFile)
        reader.Update()
        data = reader.GetOutput()

        # fill data
        # step 1, end coordinates
        nObj = int(data.GetPoints().GetNumberOfPoints() / 2)
        # print("parsing data for ", nObj, " sylinders")
        for i in range(nObj):
            s = objType()
            s.end0 = data.GetPoints().GetPoint(2 * i)
            s.end1 = data.GetPoints().GetPoint(2 * i + 1)
            objList.append(s)

        # step 2, member cell data
        numCellData = data.GetCellData().GetNumberOfArrays()
        # print("Number of CellDataArrays: ", numCellData)
        for i in range(numCellData):
            cdata = data.GetCellData().GetArray(i)
            dataName = cdata.GetName()
            # print("Parsing Cell Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName, cdata.GetTuple(j))

        # step 3, member point data
        numPointData = data.GetPointData().GetNumberOfArrays()
        # print("Number of PointDataArrays: ", numPointData)
        for i in range(numPointData):
            pdata = data.GetPointData().GetArray(i)
            dataName = pdata.GetName()
            # print("Parsing Point Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName + "0", pdata.GetTuple(2 * j))
                setattr(objList[j], dataName + "1", pdata.GetTuple(2 * j + 1))

        # output all data for debug
        # for s in objList[:10]:
        #     # print(s.end0, s.end1)
        #     attrs = vars(s)
        #     print('*************************************')
        #     print('\n'.join("%s: %s" % item for item in attrs.items()))
        #     print('*************************************')

        # print("-------------------------------------")
        return

    def parseConBlockFile(self, conBlockFile):
        self.parseFile(conBlockFile, ConBlock, self.conBlocks)
        self.conBlocks.sort(key=lambda x: (
            x.bilateral[0], x.gid0[0], x.gid1[0]))


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


def check_link(frame):
    # 4 procs, each proc has 1 chain of 20 segments so 19 bi constraints
    # block 0 : 0,1
    # block 1 : 1,2
    # block 2 : 2,3
    # .......

    biblocks = []
    for b in frame.conBlocks:
        if b.bilateral[0] > 0:
            biblocks.append(b)

    b = biblocks[0]
    chains = [[(int(b.gid0[0]), int(b.gid1[0]))]]
    biblocks.remove(biblocks[0])

    # step 1, find segments
    while biblocks:
        found = False
        for b in biblocks:
            if chains[-1][-1][1] == b.gid0[0]:
                chains[-1].append((int(b.gid0[0]), int(b.gid1[0])))
                biblocks.remove(b)
                found = True
        if not found:  # beginning of a new chain
            b = biblocks[0]
            chains.append([(int(b.gid0[0]), int(b.gid1[0]))])
            biblocks.remove(b)

    # step 2, connect found segments
    chains2 = [chains[0]]
    chains.remove(chains[0])

    while chains:
        found = False
        for c in chains:
            if chains2[-1][-1][1] == c[0][0]:
                chains2[-1].extend(c)
                chains.remove(c)
                found = True
            elif chains2[-1][0][0] == c[-1][1]:
                c.extend(chains2[-1])
                chains2[-1] = c
                chains.remove(c)
                found = True
        if not found:
            chains2.append(chains[-1])
            chains.remove(chains[-1])

    # step 3, verify chains
    assert len(chains2) == 4
    for c in chains2:
        assert len(c) == 19

    return


# get file list
# sort as numerical order
ConBlockFileList = glob.glob(
    './result/result*/ConBlock_*.pvtp')
ConBlockFileList.sort(key=getFrameNumber_lambda)

# example
frames = []
for i in range(len(ConBlockFileList)):
    # get frame
    frames.append(Frame(ConBlockFileList[i]))

for frame in frames[11:20]:
    check_link(frame)
