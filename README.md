# Integrated Sky Shader

Welcome to preintegrated sky shader with no raymarching.
That was a proof of concept that worked quite well. - to see if sky can be approximated via single formula.

Findings:
1. Yes we can make preintegrated sky
2. Mathcad cannot handle integrals so they should be replaced with more simple lerps.
3. Mie scattering is not working properly, but absolutely doable with single profile

Next steps to do:
1. Make atmospheric composition with particle profiles in script and pass single profile to shader. has to work.
2. compare with raymarching methods, it may turn out it was not worth the hassle all along.

<img width="549" alt="image" src="https://user-images.githubusercontent.com/5610313/163669812-a645b047-9e3e-4bf9-b6aa-c1683fca0557.png">
