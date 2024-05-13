import cv2
import numpy as np
from netgen.geom2d import SplineGeometry
from ngsolve import *

def createFrameGeometry(width, height, resolution=0.05):
    geo = SplineGeometry()
    geo.AddRectangle((0., 0.), (height, width))
    mesh = Mesh(geo.GenerateMesh(maxh=resolution))
    return mesh

def computeFrameGeometry(frame):
    width = frame.shape[0]
    height = frame.shape[1]
    d = np.max([width, height])
    width = width / d
    height = height / d
    return width, height
def readVideoToGridMultiDimFunction(videopath, gf, max_frames=10):
    cap = cv2.VideoCapture(videopath)
    ret = True
    count = 0
    width, height = getVideoGeometry(videopath)
    ImageSpace = gf.space
    while ret == True:
        ret, frame_new = cap.read()
        if ret == True:
            gray_img = cv2.cvtColor(frame_new, cv2.COLOR_BGR2GRAY)/1.
            func0 = VoxelCoefficient((0, 0), (height, width), gray_img, linear=True)
            func0.Compile()

            gridf = GridFunction(ImageSpace)
            gridf.Set(func0)
            gf.AddMultiDimComponent(gridf.vec)
            count += 1
            if max_frames != -1:
                if count > max_frames:
                    ret = False
    return count - 1

def getVideoGeometry(videopath):
    cap = cv2.VideoCapture(videopath)
    ret, frame = cap.read()
    return computeFrameGeometry(frame)


class NGSVideoRead():
    def __init__(self):
        self.videopath = '/home/kai/output.avi'
        self.frame_count = 10

        self.width, self.height = getVideoGeometry(self.videopath)
        self.mesh = createFrameGeometry(self.width, self.height, 0.01)

        self.ImageSpace = H1(self.mesh, order=2)
        self.FlowSpace = VectorH1(self.mesh, order=2)

        self.f = GridFunction(self.ImageSpace, multidim=0)
        self.frame_count = readVideoToGridMultiDimFunction(self.videopath, self.f, self.frame_count)

        self.gfu = GridFunction(self.ImageSpace, multidim=self.frame_count)
        for k in range(self.frame_count - 1):
            self.gfu.vecs[k].data = self.f.vecs[k]
