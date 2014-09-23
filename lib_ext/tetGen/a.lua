-- Lua script.
p=tetview:new()
p:load_mesh("C:/Users/kmmoerman/00_WORK/05_MATLAB/01_TEMP/tetGenModel_p.1.ele")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
