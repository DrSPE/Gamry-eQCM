# Garmin-eQCM
Scripts to import and resort Data from the Garmin file format into a better analysible foramt

A collection of scripts to convert the output from the Garmin eQCM setup (Garmin QCM + Potentiostat) into a file that can be plotted easily 
without having to use Garmins plotting software first.

Contains a script to sort files into their own folder.
The QCM import script splits the Garmin file format which contains the frequency data and the echem data (CV only) in different blocks and unequal spacing
and merges the same data point into the same row while interpolating missing points to make sure they are equidistant.
