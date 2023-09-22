import sys
import vtk
import glob
import re
import numpy as np
import scipy as sp
import yaml

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


def check_gamma(frame, next):
    bi_now = []
    for b in frame.conBlocks:
        if b.bilateral[0] > 0:
            bi_now.append(b)

    bi_next = []
    for b in next.conBlocks:
        if b.bilateral[0] > 0:
            bi_next.append(b)

    assert len(bi_now) == len(bi_next)
    nbi = len(bi_now)
    gamma_list = np.zeros(nbi)
    deltaKappa_list = np.zeros(nbi)
    for i in range(nbi):
        gamma_list[i] = bi_now[i].gamma[0]
        deltaKappa_list[i] = -bi_next[i].delta0[0]*linkKappa
    error = np.linalg.norm(gamma_list-deltaKappa_list) / \
        np.linalg.norm(gamma_list)
    print(error)
    if error > 2e-3:
        print("Fail: Error")

    pass


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

for i in range(len(frames)-1):
    frame = frames[i]
    next = frames[i+1]
    check_gamma(frame, next)
