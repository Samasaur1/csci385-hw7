# build system

the same as previous projects, although I have finally figured out the issue that TypeScript was complaining about, so I don't get errors anymore

# the project

This felt pretty straightforward --- just follow the instructions and implement the given math. The part that was most challenging for me is that until I implemented Verlet integration, nothing in my simulation worked: the sheet would just spaz out and shoot offscreen. And since I did that last, I thought that I had other problems in my code until I implemented it and suddenly everything started working.

For some reason, my simulation is super slow. I have adjusted `gTimeStep` to make it slightly faster, but if it's too fast it becomes really unstable. I also adjusted `gStiffness` a little.

## possible extension

I noticed when trying it out that the wind doesn't *really* emanate from the propeller of the plane. Maybe it would be interesting to only have masses be affected by the wind if they are below a certain height? or affected by the wind at different strengths based on height? I imagine this would probably produce a fluttering effect at the end of the sheet, as opposed to the current behavior of eventually reaching a steady-state position.
