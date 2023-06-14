# evDisp

evDisp

## Data

Describe any data formats evDisp creates, destroys, changes, analyzes, or its usage.

## Configuration

Describe any configuration variables for evDisp.

```
param1 value1
param2 value2
```

## Functions

### Initialise
- Instances the 2D graphs/histograms that are used to plot the top, bottom and barrel of the detector.
- Instances the histograms sued to plot the hits against time and the hits against charge.
- Instances the canvas and divides it into sections so that the top and bottom of the barrel can be plotted alongside each other with the barrel below and the charge/time histograms alongside each other below the barrel.
- Sets up the colour palette to use on the graphs that is used to plot the PMT charge distribution on the event display. The colour pallete emulates the following p.e ranges:
  - > 26.7     rgb(196, 032, 033)
  - 23.3-26.7  rgb(219, 109, 040)
  - 20.2-23.3  rgb(210, 148, 046)
  - 17.3-20.2  rgb(220, 185, 055)
  - 14.7-17.3  rgb(223, 187, 061)
  - 12.2-14.7  rgb(236, 240, 072)
  - 10.0-12.2  rgb(245, 252, 079)
  - 08.0-10.0   rgb(206, 226, 067)
  - 06.2-08.0    rgb(165, 186, 052)
  - 04.7-06.2    rgb(146, 178, 052)
  - 03.3-04.7    rgb(100, 140, 050)
  - 02.2-03.3    rgb(094, 137, 117)
  - 01.3-02.2    rgb(088, 107, 215)
  - 00.7-01.3    rgb(029, 032, 205)
  - 00.2-00.7    rgb(084, 043, 209)
  - < 0.2      rgb(121, 053, 207)
- Note that ROOT does not plot over the full colour range for some reason so the palette values have been adjusted down to fix this.

### Execute
- Pull hit information from the relevant location as defined in the config file.
- The common block sktqz_.
- The branch TQREAL.
- A combination of the common blocks skchnl_, skq_ and skt_.
- Note that event information (run number, subrun number, event number and trigger ID) is pulled from the either the skhead_ common block or the Header branch.
