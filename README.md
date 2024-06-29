# dark_matter_demo
Numerical simulation showing that the observed galactic mass is not enough to justify the velocities of stars for from the galactic center.

This script calculates the gravitational force and velocity profile for a galaxy with a given mass density distribution. This demonstration uses
arbitrary values and a numerical approach to demonstrate how our current understanding of gravity and mass does not align with our observations.
In reality, there is a theory that dark matter changes the mass distribution of galaxies to allow for farther objects to travel faster. 

The model generates a galaxy of 'stars' by placing point-masses along many concentric circles of increasing diameter. The galaxy is flat, so
there is no out-of-plane forces being modeled. Along the x-axis, where y=0, the force on that particular point is towards the center of the
galaxy. The magnitude will depend on the mass distribution. The force in the y-direction is zero, because the mass is symmetric about the x-axis.
Then, using classical mechanics, the velocity required to keep a star in a circular orbit is calculated. This is what we should expect during
observations.

The user can tweak the parameters in the first code block below to simulate different mass density distribution equations or use
non-arbitrary values to get real numbers.

