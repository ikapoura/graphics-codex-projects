# Graphics Codex projects

This is my attempt to implement the projects by Morgan Mcguire in https://graphicscodex.com/projects/projects/index.html. It is a good refresher for the basics and gives concrete goals that is excellent when distracted by many things.

Since I am actively working on it, new information will surface as I push through them. My notes and a simple review of each project is written in [my blog](https://iliaskapouranis.com/2021/04/24/studying-with-the-graphics-codex/).

## 1. Cubes
The first project's goal is to create a scene using only cubes. I implemented a way to read an image and create a cube grid for each pixel, where each cube's height is determined by the corresponding pixel's luminance. Simple pixels to cubes.

|<img width="1024" height="300" alt="Horizontal grid" src="https://user-images.githubusercontent.com/40468844/115958746-512f1800-a511-11eb-89cc-8dd064a8235f.jpg"> Horizontal view | <img width="1024" height="300" alt="Vertical grid" src="https://user-images.githubusercontent.com/40468844/115958748-52604500-a511-11eb-8e81-2ed0289c0829.jpg"> Vertical view |
|:-:|:-:|

## 2. Meshes
For the second project, we move to something more ambitious where we have to procedurally generate a wine glass. This can be done however we want, but I decided to utilize images again for the generation. The user needs to provide a grayscale image which contains the half of the cross-section of the model they wish to generate. Then, using a "Quality" slider, they can change how much geometry is generated and how faithful the final mesh is to the input image. The "Quality" slider is using a distance and an angle criterion in order to reject triangle rings. The lower the "Quality" the less geometry is generated, although one can notice that the low "Quality" meshes are cleaner.!

|<img width="200" height="300" alt="Input grayscale image" src="https://user-images.githubusercontent.com/40468844/115959165-8a688780-a513-11eb-9107-597f0f01d017.png"> <p> Input grayscale image</p> | <img width="256" height="300" alt="Low quality" src="https://user-images.githubusercontent.com/40468844/115959170-92c0c280-a513-11eb-9e3b-beaf76f801c4.png"> <p>Low quality without distance and angle criteria</p> | <img width="220" height="300" alt="Low quality with distance and angle criteria" src="https://user-images.githubusercontent.com/40468844/115959175-9d7b5780-a513-11eb-9101-a82b7b19d8f4.png"> <p>Low quality with distance and angle criteria</p> |
|:-:|:-:|:-:|

After settling to a solution that is good enough, I had to make a scene to present the results in an appealing way. For this reason, I created a more refined wine glass and a plate model with the same technique and arranged them in a bar setting. Below are the input image and resulting mesh for the wine glass and the plates and a render of them in a bar skybox.

|<img width="100" height="300" alt="Wine glass grayscale image input" src="https://user-images.githubusercontent.com/40468844/115959519-49717280-a515-11eb-815a-85b627757d6f.png"> <p>Wine glass grayscale image input</p> | <img width="210" height="300" alt="Wine glass mesh output" src="https://user-images.githubusercontent.com/40468844/115959532-58582500-a515-11eb-90b7-e8c7a50562c3.png"> <p>Wine glass mesh output</p> |
|:-:|:-:|
|<img width="600" height="300" alt="Plate grayscale image input" src="https://user-images.githubusercontent.com/40468844/115959524-51311700-a515-11eb-9877-72e2134465ad.png"> <p>Plate grayscale image input</p> | <img width="640" height="300" alt="Plate mesh output" src="https://user-images.githubusercontent.com/40468844/115959533-5a21e880-a515-11eb-812a-2caa8fb90a4d.png"> <p>Plate mesh output</p> |

|<img width="1024" height="562" alt="Rendering of the wine glass and plate in G3D" src="https://user-images.githubusercontent.com/40468844/115959564-74f45d00-a515-11eb-8f4c-7ec32a47a794.png"> <p>Rendering of the wine glass and plate in G3D</p> |
|:-:|

