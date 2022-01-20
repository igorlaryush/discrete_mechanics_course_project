import numpy as np
import matplotlib.pyplot as plt

from celluloid import Camera
fig, ax = plt.subplots(figsize=(10,10))
camera = Camera(fig)

from Consts import *
from Point import point
from Spring import  spring


points = [None] * POINTS_NUMBER
springs = [None] * NUMS
wall = [None] * WALL_NUMS


def CreateBall(C):
    for i in range(POINTS_NUMBER):
        points[i] = point(C.x + RADIUS * np.cos(i * 2 * np.pi / POINTS_NUMBER), C.y + RADIUS * np.sin(i * 2 * np.pi / POINTS_NUMBER), v0x, v0y, 0, 0)
    for i in range(POINTS_NUMBER - 1):
        connect_spring(i, i, i + 1)
    connect_spring(POINTS_NUMBER - 1, POINTS_NUMBER - 1, 0)


def DrawBall(iter):
    x_w = np.array([wall[i].x  for i in range(WALL_NUMS)])
    y_w = np.array([wall[i].y  for i in range(WALL_NUMS)])
    x = [points[i].x for i in range(POINTS_NUMBER)]
    x.append(points[0].x)
    x = np.array(x)
    y = [points[i].y for i in range(POINTS_NUMBER)]
    y.append(points[0].y)
    y = np.array(y)
    ax.plot(x, y, "r")
    camera.snap()


def CalcForces():
    for i in range(POINTS_NUMBER):
        points[i].fx = 0.0
        points[i].fy = 0.0
    # spring forces
    for i in range(POINTS_NUMBER):
        x1 = points[springs[i].i].x
        y1 = points[springs[i].i].y
        x2 = points[springs[i].j].x
        y2 = points[springs[i].j].y
        r12d = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        if r12d != 0:
            vx12 = points[springs[i].i].vx - points[springs[i].j].vx
            vy12 = points[springs[i].i].vy - points[springs[i].j].vy

            f = (r12d - springs[i].length) * KS + (vx12 * (x1 - x2) + vy12 * (y1 - y2)) * KD / r12d

            Fx = ((x1 - x2) / r12d) * f
            Fy = ((y1 - y2) / r12d) * f

            points[springs[i].i].fx -= Fx
            points[springs[i].i].fy -= Fy

            points[springs[i].j].fx += Fx
            points[springs[i].j].fy += Fy
        springs[i].nx = -(y1 - y2) / r12d
        springs[i].ny = (x1 - x2) / r12d
    volume = 0.0
    for i in range(POINTS_NUMBER - 1):
        x1 = points[springs[i].i].x
        y1 = points[springs[i].i].y
        x2 = points[springs[i].j].x
        y2 = points[springs[i].j].y
        r12d = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        volume += 0.5 * np.abs(x1 - x2) * np.abs(springs[i].nx) * (r12d)
    # pressure
    for i in range(POINTS_NUMBER):
        x1 = points[springs[i].i].x
        y1 = points[springs[i].i].y
        x2 = points[springs[i].j].x
        y2 = points[springs[i].j].y
        r12d = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        pressurev = r12d * Pressure * (1.0 / volume)
        points[springs[i].i].fx += springs[i].nx * pressurev
        points[springs[i].i].fy += springs[i].ny * pressurev
        points[springs[i].j].fx += springs[i].nx * pressurev
        points[springs[i].j].fy += springs[i].ny * pressurev
    # wall
    for i in range(POINTS_NUMBER):
        x = points[i].x
        y = points[i].y
        f_w_x = 0.0
        f_w_y = 0.0
        r = np.abs(x)
        gradLD = grad(r)
        points[i].fx += -gradLD


def integrate():
    vt_2_x = []
    vt_2_y = []
    for i in range(POINTS_NUMBER):
        points[i].x += points[i].vx * DT + 0.5 * (points[i].fx / MASS )* DT * DT
        vt_2_x.append(points[i].vx + 0.5 * (points[i].fx / MASS )* DT)
        points[i].y += points[i].vy * DT + 0.5 * (points[i].fy / MASS )* DT * DT
        vt_2_y.append(points[i].vy + 0.5 * (points[i].fy / MASS )* DT)
    CalcForces()
    for i in range(POINTS_NUMBER):
        points[i].vx = vt_2_x[i] +  0.5 * (points[i].fx / MASS )* DT
        points[i].vy = vt_2_y[i] +  0.5 * (points[i].fy / MASS )* DT


def grad(r):
    e = 10**(-10)
    a = 1
    return 24 * a**6 * e * (-2 * a**6 + r**6) / r**13

def Wall():
    h = WALL_HEIGHT / WALL_NUMS
    for i in range(WALL_NUMS):
        wall[i] = point(0, i * h, 0, 0, 0, 0)


def connect_spring(pi, i, j):
    springs[pi] = spring(0, 0, 0, 0, 0)
    springs[pi].i = i;
    springs[pi].j = j;
    springs[pi].length = np.sqrt(
       (points[ i ].x - points[ j ].x) *
       (points[ i ].x - points[ j ].x) +
       (points[ i ].y - points[ j ].y) *
       (points[ i ].y - points[ j ].y)
      )


c = point(5.0, 3.0, 0, 0, 0, 0)
CreateBall(c)
Wall()
CalcForces()
for i in range(iters):
    integrate()
    if Pressure < FINAL_PRESSURE:
        Pressure += FINAL_PRESSURE / 50.0
    if i % 50 == 0:
        DrawBall(i)
anim = camera.animate()