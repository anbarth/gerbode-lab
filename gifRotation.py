import imageio
import os

images = []

#for i in range(50):
#    images.append(imageio.imread('pic1.png'))
#for i in range(50):
#    images.append(imageio.imread('pic2.png'))

directory = 'hell'
for entry in os.scandir(directory):
    if (entry.path.endswith(".jpg")
            or entry.path.endswith(".png")) and entry.is_file():
        #images.append(entry)
        images.append(imageio.imread(entry.path))
        #print(entry.path)

# make the gif!
imageio.mimsave('gifs/waaa.gif', images)