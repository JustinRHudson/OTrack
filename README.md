OLR Track (OTrack) began as a way to track convective Outgoing Longwave
Radiation (OLR) anomalies related to the Madden-Julian Oscillation (MJO).
The original version (preserved in this repo as OLR_Tracker.py) was used
in my 2023 paper "The Role of Surface Fluxes in MJO Propagation Through
The Maritime Continent". 

Since then I have been working on improving the core functionality and
efficiency of this code. The current version is Blob_Tracker.py and is a
work in progress. The original OLR_Tracker is not limited to just tracking
OLR anomalies (see cheese.gif and dvd.gif under OLD_EXAMPLES) but has some
hard coded features specific to the peculiarities of my MJO tracking.

The original OLR_Tracker mainly uses base Python (with a little bit of numpy).
Blob_Tracker.py continues in this vein of being package minimal, although as I
encounter problems that I believe have been best solved elsewhere I plan on
implementing other publically available packages.

Utilizing the tracking packages in this repo is simple. Simply import the .py
file from its location on your machine and uses the Track function on a binary
map you have created. The code assumes a format of (time,lat,lon) The code will
then return a dictionary of Blob Objects (I like to call them Blobjects) which
contain position, velocity, and temporal data for all Blobs tracked in your
binary map.