import numpy as np
import cv2

# initialize water image
height = 300
width = 300
water_depth = np.zeros((height, width), dtype=float)

# initialize video writer
fourcc = cv2.VideoWriter_fourcc('M','J','P','G')
fps = 30
video_filename = 'output.avi'
out = cv2.VideoWriter(video_filename, fourcc, fps, (width, height))

# new frame after each addition of water
for i in range(width - 2*50):
    #random_locations = np.random.random_integers(200,450, size=(200, 2))
    #for item in random_locations:
    #water_depth[item[0], item[1]] += 0.1
    frame = np.zeros((height, width))
    frame = cv2.circle(frame, (int(50) + i, int(height/2)), 50, (1., 0., 0.), -1)
    #add this array to the video
    gray = cv2.normalize(frame, None, 255, 0, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    gray_3c = cv2.merge([gray, gray, gray])
    out.write(gray_3c)

# close out the video writer
out.release()