# removalBias

The files in the repository were used to calculate the removal bias estimates for the [NetAssessApp](https://github.com/LADCO/NetAssessApp).

# Background
The US EPA requires state agencies to assess their air pollution monitoring networks every 5 years (40 CFR part 58.10(d)). The [NetAssessApp](https://github.com/LADCO/NetAssessApp) was developed by [LADCO](http://www.ladco.org/) and state enironmental agencies in EPA Region 5 as a resource for doing network assessment analyses. The online tool was designed as an attempt to update the network assessment tools provided by the EPA for the 2010 network assessment ([see EPA's site for network assessment tools and guidance](http://www.epa.gov/ttnamti1/network-assessment.html)).

# Bias Estimates
The removal bias tool is meant to aid in determining redundant sites. The bias estimation uses the nearest neighbors to each site to estimate the concentration at the location of the site if the site had never existed.  This is done using the Voronoi Neighborhood Averaging algorithm with inverse distance squared weighting. The squared distance allows for higher weighting on concentrations at sites located closer to the site being examined.  The bias was calculated for each day at each site by taking the difference between the predicted value from the interpolation and the measured concentration. A positive average bias would mean that if the site being examined was removed, the neighboring sites would indicate that the estimated concentration would be larger than the measured concentration.  Likewise, a negative average bias would suggest that the estimated concentration at the location of the site is smaller than the actual measured concentration. 
