OLR Track (OTrack) began as a way to track convective Outgoing Longwave
Radiation (OLR) anomalies related to the Madden-Julian Oscillation (MJO).
The original version (preserved in this repo as OLR_Tracker.py) was used
in my 2023 paper "The Role of Surface Fluxes in MJO Propagation Through
The Maritime Continent". 

Since then I have been working on improving the core functionality and
efficiency of this code. The current version is Blob_Tracker.py and is a
work in progress (The current version on the repo is not fully updated
and this script should be considered experimental). The original OLR_Tracker
has some hard coded features specific to the peculiarities of my MJO tracking.

The original OLR_Tracker mainly uses base Python (with a little bit of numpy).
Blob_Tracker.py continues in this vein of being package minimal, although as I
encounter problems that I believe have been best solved elsewhere I plan on
implementing other publically available packages.

Utilizing the tracking packages in this repo is simple. Simply import the .py
file from its location on your machine and uses the Track function on a binary
map you have created. The code assumes a format of (time,lat,lon) The code will
then return a dictionary of tracked objects which contain position, velocity, 
and temporal data for all Blobs tracked from your binary map.