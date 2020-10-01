
The code
[get_scenarios_similar_to_historical_events.R](get_scenarios_similar_to_historical_events.R)
shows how to download and plot the initial conditions for PTHA18 scenarios that
had top-3 goodness of fit when compared with historical data at DART buoys. In the
case of VAUS scenarios we include more than 3 if there are "double-ups" (i.e. identical
scenarios, which can easily happen with uniform-slip).

It uses a hard-coded record of the 'best scenarios' that was created with
[find_desired_event_rows.R](find_desired_event_rows.R), and is stored in (for
heterogeneous-slip) [best_fitting_HS.R](best_fitting_HS.R) and (for
variable-area-uniform-slip) [best_fitting_VAUS.R](best_fitting_VAUS.R). 

There is no guarentee that these scenarios will match the tsunami as observed
elsewhere -- sometimes they will, but not always. In principle one should
expect to get better accuracy using an inverted source that considers data
at your specific site of interest.

*2020/10/02: I noticed that in some VAUS cases we are finding "near
double-ups" with slip differing by 0.001m (see the VAUS cases for the puysegur2009
scenario). Currently such VAUS "near double-ups" are treated as distinct,
although they probably shouldn't be. I was surprised that such "near double-up"
scenarios exist - however on reflection, the code to make the vaus-scenarios
sets the slip based on a corresponding heterogeneous-slip scenario by summing
its (slip x area) and dividing by the total area. With different
heterogeneous-slip, there is potential for round-off level changes in the sum
of (slip x area). Our post-hoc tests confirmed that the integrated (slip x area
x rigidity) of all PTHA18 scenarios was within 1/1000 of the target value,
so any such differences will necessarily be very small and unimportant in regards
of the earthquake magnitude.*
