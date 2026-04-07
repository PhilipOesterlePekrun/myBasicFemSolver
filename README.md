This is a personal project of mine, a sort of proof of concept of what I know about FEM, from the ground up, with no external libraries except for [SFML](https://github.com/SFML/SFML) for visualization (and my own libraries). I tend to move most of the general methods (e.g. numerical integration, linear solvers) into [my utilities library](https://github.com/PhilipOesterlePekrun/myUtils), so myBasicFemSolver makes extensive use of it.

The focus is FEM for solid mechanics problems, but I am leaving room for other physics and even other discretization methods or domains of computational mechanics (e.g. rigid body, free surface flow), as well as things like meshing and boundary condition definition. Currently, linear FEM for static and dynamic (using Newmark-beta method) 2D solid mechanics works very well. Here are some random pictures and gifs:

Youngs = Density = Horizontal length = 1.0; Poissons = 0.4:
![recording_compressed](https://github.com/user-attachments/assets/f4020046-222c-4207-a0cc-16fdaf69b519)

Material properties =~ Aluminium; Horizontal length = 1.0:
![recording2_compressed](https://github.com/user-attachments/assets/e72b2ee1-eb1a-400e-b04b-bc77c1ecd3da)

Static:
<img width="1600" height="1032" alt="Screenshot_20260309_184635" src="https://github.com/user-attachments/assets/7c4eba65-d8fb-43ab-92c3-c28ced659e78" />
<img width="1600" height="1026" alt="Screenshot_20260304_233308" src="https://github.com/user-attachments/assets/fe1054f1-3d77-447c-b527-1d9002efa31f" />
<img width="1601" height="1034" alt="Screenshot_20260309_192147" src="https://github.com/user-attachments/assets/49182b6e-2562-48ad-b14a-746cc63c410f" />
<img width="1605" height="1028" alt="Screenshot_20260309_192127" src="https://github.com/user-attachments/assets/3649efd0-224c-4fee-a0a9-b509655877df" />
