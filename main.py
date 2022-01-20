from Functions import *
import matplotlib

fig, ax = plt.subplots(figsize=(20,20))
camera = Camera(fig)

writergif = matplotlib.animation.PillowWriter(fps=15)
anim.save('filename.gif', writer=writergif)