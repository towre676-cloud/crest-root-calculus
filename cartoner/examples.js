globalThis.trialExamples = {
  "opener": {
    "setup": {
      "schema": 1,
      "setup": "opener-r1",
      "name": "Carton opening and transfer",
      "source": "synthetic",
      "channels": [
        {
          "id": "opening",
          "label": "Opening angle",
          "unit": "deg",
          "sense": "min",
          "limit": 86,
          "error": 0.5,
          "rate": 1.2,
          "window": [
            160,
            180
          ]
        },
        {
          "id": "vacuum",
          "label": "Vacuum magnitude",
          "unit": "kPa",
          "sense": "min",
          "limit": 38,
          "error": 1,
          "rate": 0.6,
          "window": [
            120,
            170
          ]
        }
      ],
      "events": [
        {
          "id": "captured",
          "errorDeg": 0.8
        },
        {
          "id": "released",
          "errorDeg": 0.6
        }
      ],
      "handoffs": [
        {
          "label": "Capture before vacuum release",
          "first": "captured",
          "second": "released",
          "minDeg": 4
        }
      ],
      "loops": [
        {
          "name": "Reference carrier crank",
          "groundDeg": 0,
          "ground": [
            93.55,
            94.05
          ],
          "crank": [
            46.35,
            46.85
          ],
          "coupler": [
            83.85,
            84.35
          ],
          "rocker": [
            111.05,
            111.55
          ]
        }
      ]
    },
    "csv": "setup,run,cycle,phase_deg,channel,value\nopener-r1,synthetic-bench,1,0,opening,0.0\nopener-r1,synthetic-bench,1,0,vacuum,8.0\nopener-r1,synthetic-bench,1,10,opening,0.0\nopener-r1,synthetic-bench,1,10,vacuum,12.2105\nopener-r1,synthetic-bench,1,20,opening,0.0\nopener-r1,synthetic-bench,1,20,vacuum,16.4211\nopener-r1,synthetic-bench,1,30,opening,0.0\nopener-r1,synthetic-bench,1,30,vacuum,20.6316\nopener-r1,synthetic-bench,1,40,opening,0.0\nopener-r1,synthetic-bench,1,40,vacuum,24.8421\nopener-r1,synthetic-bench,1,50,opening,8.0455\nopener-r1,synthetic-bench,1,50,vacuum,29.0526\nopener-r1,synthetic-bench,1,60,opening,16.0909\nopener-r1,synthetic-bench,1,60,vacuum,33.2632\nopener-r1,synthetic-bench,1,70,opening,24.1364\nopener-r1,synthetic-bench,1,70,vacuum,37.4737\nopener-r1,synthetic-bench,1,80,opening,32.1818\nopener-r1,synthetic-bench,1,80,vacuum,41.6842\nopener-r1,synthetic-bench,1,90,opening,40.2273\nopener-r1,synthetic-bench,1,90,vacuum,45.8947\nopener-r1,synthetic-bench,1,100,opening,48.2727\nopener-r1,synthetic-bench,1,100,vacuum,48.0\nopener-r1,synthetic-bench,1,110,opening,56.3182\nopener-r1,synthetic-bench,1,110,vacuum,48.0\nopener-r1,synthetic-bench,1,120,opening,64.3636\nopener-r1,synthetic-bench,1,120,vacuum,48.0\nopener-r1,synthetic-bench,1,130,opening,72.4091\nopener-r1,synthetic-bench,1,130,vacuum,48.0\nopener-r1,synthetic-bench,1,140,opening,80.4545\nopener-r1,synthetic-bench,1,140,vacuum,48.0\nopener-r1,synthetic-bench,1,150,opening,88.5\nopener-r1,synthetic-bench,1,150,vacuum,48.0\nopener-r1,synthetic-bench,1,160,opening,88.5\nopener-r1,synthetic-bench,1,160,vacuum,48.0\nopener-r1,synthetic-bench,1,170,opening,88.5\nopener-r1,synthetic-bench,1,170,vacuum,48.0\nopener-r1,synthetic-bench,1,180,opening,88.5\nopener-r1,synthetic-bench,1,180,vacuum,48.0\nopener-r1,synthetic-bench,1,190,opening,88.5\nopener-r1,synthetic-bench,1,190,vacuum,48.0\nopener-r1,synthetic-bench,1,200,opening,81.125\nopener-r1,synthetic-bench,1,200,vacuum,44.3636\nopener-r1,synthetic-bench,1,210,opening,73.75\nopener-r1,synthetic-bench,1,210,vacuum,40.7273\nopener-r1,synthetic-bench,1,220,opening,66.375\nopener-r1,synthetic-bench,1,220,vacuum,37.0909\nopener-r1,synthetic-bench,1,230,opening,59.0\nopener-r1,synthetic-bench,1,230,vacuum,33.4545\nopener-r1,synthetic-bench,1,240,opening,51.625\nopener-r1,synthetic-bench,1,240,vacuum,29.8182\nopener-r1,synthetic-bench,1,250,opening,44.25\nopener-r1,synthetic-bench,1,250,vacuum,26.1818\nopener-r1,synthetic-bench,1,260,opening,36.875\nopener-r1,synthetic-bench,1,260,vacuum,22.5455\nopener-r1,synthetic-bench,1,270,opening,29.5\nopener-r1,synthetic-bench,1,270,vacuum,18.9091\nopener-r1,synthetic-bench,1,280,opening,22.125\nopener-r1,synthetic-bench,1,280,vacuum,15.2727\nopener-r1,synthetic-bench,1,290,opening,14.75\nopener-r1,synthetic-bench,1,290,vacuum,11.6364\nopener-r1,synthetic-bench,1,300,opening,7.375\nopener-r1,synthetic-bench,1,300,vacuum,8.0\nopener-r1,synthetic-bench,1,310,opening,0.0\nopener-r1,synthetic-bench,1,310,vacuum,8.0\nopener-r1,synthetic-bench,1,320,opening,0.0\nopener-r1,synthetic-bench,1,320,vacuum,8.0\nopener-r1,synthetic-bench,1,330,opening,0.0\nopener-r1,synthetic-bench,1,330,vacuum,8.0\nopener-r1,synthetic-bench,1,340,opening,0.0\nopener-r1,synthetic-bench,1,340,vacuum,8.0\nopener-r1,synthetic-bench,1,350,opening,0.0\nopener-r1,synthetic-bench,1,350,vacuum,8.0\nopener-r1,synthetic-bench,1,360,opening,0.0\nopener-r1,synthetic-bench,1,360,vacuum,8.0\nopener-r1,synthetic-bench,1,172,captured,\nopener-r1,synthetic-bench,1,181,released,\nopener-r1,synthetic-bench,2,0,opening,0.0\nopener-r1,synthetic-bench,2,0,vacuum,8.0\nopener-r1,synthetic-bench,2,10,opening,0.0\nopener-r1,synthetic-bench,2,10,vacuum,12.2105\nopener-r1,synthetic-bench,2,20,opening,0.0\nopener-r1,synthetic-bench,2,20,vacuum,16.4211\nopener-r1,synthetic-bench,2,30,opening,0.0\nopener-r1,synthetic-bench,2,30,vacuum,20.6316\nopener-r1,synthetic-bench,2,40,opening,0.0\nopener-r1,synthetic-bench,2,40,vacuum,24.8421\nopener-r1,synthetic-bench,2,50,opening,8.4545\nopener-r1,synthetic-bench,2,50,vacuum,29.0526\nopener-r1,synthetic-bench,2,60,opening,16.9091\nopener-r1,synthetic-bench,2,60,vacuum,33.2632\nopener-r1,synthetic-bench,2,70,opening,25.3636\nopener-r1,synthetic-bench,2,70,vacuum,37.4737\nopener-r1,synthetic-bench,2,80,opening,33.8182\nopener-r1,synthetic-bench,2,80,vacuum,41.6842\nopener-r1,synthetic-bench,2,90,opening,42.2727\nopener-r1,synthetic-bench,2,90,vacuum,45.8947\nopener-r1,synthetic-bench,2,100,opening,50.7273\nopener-r1,synthetic-bench,2,100,vacuum,48.0\nopener-r1,synthetic-bench,2,110,opening,59.1818\nopener-r1,synthetic-bench,2,110,vacuum,48.0\nopener-r1,synthetic-bench,2,120,opening,67.6364\nopener-r1,synthetic-bench,2,120,vacuum,48.0\nopener-r1,synthetic-bench,2,130,opening,76.0909\nopener-r1,synthetic-bench,2,130,vacuum,48.0\nopener-r1,synthetic-bench,2,140,opening,84.5455\nopener-r1,synthetic-bench,2,140,vacuum,48.0\nopener-r1,synthetic-bench,2,150,opening,93.0\nopener-r1,synthetic-bench,2,150,vacuum,48.0\nopener-r1,synthetic-bench,2,160,opening,93.0\nopener-r1,synthetic-bench,2,160,vacuum,48.0\nopener-r1,synthetic-bench,2,170,opening,93.0\nopener-r1,synthetic-bench,2,170,vacuum,48.0\nopener-r1,synthetic-bench,2,180,opening,93.0\nopener-r1,synthetic-bench,2,180,vacuum,48.0\nopener-r1,synthetic-bench,2,190,opening,93.0\nopener-r1,synthetic-bench,2,190,vacuum,48.0\nopener-r1,synthetic-bench,2,200,opening,85.25\nopener-r1,synthetic-bench,2,200,vacuum,44.3636\nopener-r1,synthetic-bench,2,210,opening,77.5\nopener-r1,synthetic-bench,2,210,vacuum,40.7273\nopener-r1,synthetic-bench,2,220,opening,69.75\nopener-r1,synthetic-bench,2,220,vacuum,37.0909\nopener-r1,synthetic-bench,2,230,opening,62.0\nopener-r1,synthetic-bench,2,230,vacuum,33.4545\nopener-r1,synthetic-bench,2,240,opening,54.25\nopener-r1,synthetic-bench,2,240,vacuum,29.8182\nopener-r1,synthetic-bench,2,250,opening,46.5\nopener-r1,synthetic-bench,2,250,vacuum,26.1818\nopener-r1,synthetic-bench,2,260,opening,38.75\nopener-r1,synthetic-bench,2,260,vacuum,22.5455\nopener-r1,synthetic-bench,2,270,opening,31.0\nopener-r1,synthetic-bench,2,270,vacuum,18.9091\nopener-r1,synthetic-bench,2,280,opening,23.25\nopener-r1,synthetic-bench,2,280,vacuum,15.2727\nopener-r1,synthetic-bench,2,290,opening,15.5\nopener-r1,synthetic-bench,2,290,vacuum,11.6364\nopener-r1,synthetic-bench,2,300,opening,7.75\nopener-r1,synthetic-bench,2,300,vacuum,8.0\nopener-r1,synthetic-bench,2,310,opening,0.0\nopener-r1,synthetic-bench,2,310,vacuum,8.0\nopener-r1,synthetic-bench,2,320,opening,0.0\nopener-r1,synthetic-bench,2,320,vacuum,8.0\nopener-r1,synthetic-bench,2,330,opening,0.0\nopener-r1,synthetic-bench,2,330,vacuum,8.0\nopener-r1,synthetic-bench,2,340,opening,0.0\nopener-r1,synthetic-bench,2,340,vacuum,8.0\nopener-r1,synthetic-bench,2,350,opening,0.0\nopener-r1,synthetic-bench,2,350,vacuum,8.0\nopener-r1,synthetic-bench,2,360,opening,0.0\nopener-r1,synthetic-bench,2,360,vacuum,8.0\nopener-r1,synthetic-bench,2,172,captured,\nopener-r1,synthetic-bench,2,181,released,\n"
  },
  "folder": {
    "setup": {
      "schema": 1,
      "setup": "folder-r1",
      "name": "Front and rear flap folding",
      "source": "synthetic",
      "channels": [
        {
          "id": "flap",
          "label": "Rear flap angle",
          "unit": "deg",
          "sense": "min",
          "limit": 82,
          "error": 0.75,
          "rate": 0.8,
          "window": [
            240,
            260
          ]
        },
        {
          "id": "gap",
          "label": "Return clearance",
          "unit": "mm",
          "sense": "min",
          "limit": 1,
          "error": 0.1,
          "rate": 0.18,
          "window": [
            290,
            350
          ]
        }
      ],
      "events": [
        {
          "id": "rear_folded",
          "errorDeg": 1
        },
        {
          "id": "tool_clear",
          "errorDeg": 0.6
        }
      ],
      "handoffs": [
        {
          "label": "Rear fold before withdrawal",
          "first": "rear_folded",
          "second": "tool_clear",
          "minDeg": 6
        }
      ],
      "loops": [
        {
          "name": "Illustrative folding actuator",
          "groundDeg": 0,
          "ground": [
            93.55,
            94.05
          ],
          "crank": [
            46.35,
            46.85
          ],
          "coupler": [
            83.85,
            84.35
          ],
          "rocker": [
            111.05,
            111.55
          ],
          "strokeCos": [
            0,
            0.8
          ]
        }
      ]
    },
    "csv": "setup,run,cycle,phase_deg,channel,value\nfolder-r1,synthetic-bench,1,0,flap,0.0\nfolder-r1,synthetic-bench,1,0,gap,9.0\nfolder-r1,synthetic-bench,1,10,flap,0.0\nfolder-r1,synthetic-bench,1,10,gap,9.0\nfolder-r1,synthetic-bench,1,20,flap,0.0\nfolder-r1,synthetic-bench,1,20,gap,9.0\nfolder-r1,synthetic-bench,1,30,flap,0.0\nfolder-r1,synthetic-bench,1,30,gap,9.0\nfolder-r1,synthetic-bench,1,40,flap,0.0\nfolder-r1,synthetic-bench,1,40,gap,9.0\nfolder-r1,synthetic-bench,1,50,flap,0.0\nfolder-r1,synthetic-bench,1,50,gap,9.0\nfolder-r1,synthetic-bench,1,60,flap,0.0\nfolder-r1,synthetic-bench,1,60,gap,9.0\nfolder-r1,synthetic-bench,1,70,flap,0.0\nfolder-r1,synthetic-bench,1,70,gap,9.0\nfolder-r1,synthetic-bench,1,80,flap,0.0\nfolder-r1,synthetic-bench,1,80,gap,9.0\nfolder-r1,synthetic-bench,1,90,flap,0.0\nfolder-r1,synthetic-bench,1,90,gap,9.0\nfolder-r1,synthetic-bench,1,100,flap,0.0\nfolder-r1,synthetic-bench,1,100,gap,9.0\nfolder-r1,synthetic-bench,1,110,flap,0.0\nfolder-r1,synthetic-bench,1,110,gap,9.0\nfolder-r1,synthetic-bench,1,120,flap,0.0\nfolder-r1,synthetic-bench,1,120,gap,9.0\nfolder-r1,synthetic-bench,1,130,flap,7.4167\nfolder-r1,synthetic-bench,1,130,gap,9.0\nfolder-r1,synthetic-bench,1,140,flap,14.8333\nfolder-r1,synthetic-bench,1,140,gap,9.0\nfolder-r1,synthetic-bench,1,150,flap,22.25\nfolder-r1,synthetic-bench,1,150,gap,9.0\nfolder-r1,synthetic-bench,1,160,flap,29.6667\nfolder-r1,synthetic-bench,1,160,gap,9.0\nfolder-r1,synthetic-bench,1,170,flap,37.0833\nfolder-r1,synthetic-bench,1,170,gap,9.0\nfolder-r1,synthetic-bench,1,180,flap,44.5\nfolder-r1,synthetic-bench,1,180,gap,9.0\nfolder-r1,synthetic-bench,1,190,flap,51.9167\nfolder-r1,synthetic-bench,1,190,gap,9.0\nfolder-r1,synthetic-bench,1,200,flap,59.3333\nfolder-r1,synthetic-bench,1,200,gap,9.0\nfolder-r1,synthetic-bench,1,210,flap,66.75\nfolder-r1,synthetic-bench,1,210,gap,9.0\nfolder-r1,synthetic-bench,1,220,flap,74.1667\nfolder-r1,synthetic-bench,1,220,gap,9.0\nfolder-r1,synthetic-bench,1,230,flap,81.5833\nfolder-r1,synthetic-bench,1,230,gap,9.0\nfolder-r1,synthetic-bench,1,240,flap,89.0\nfolder-r1,synthetic-bench,1,240,gap,9.0\nfolder-r1,synthetic-bench,1,250,flap,89.0\nfolder-r1,synthetic-bench,1,250,gap,8.3577\nfolder-r1,synthetic-bench,1,260,flap,89.0\nfolder-r1,synthetic-bench,1,260,gap,7.0731\nfolder-r1,synthetic-bench,1,270,flap,81.3611\nfolder-r1,synthetic-bench,1,270,gap,5.7885\nfolder-r1,synthetic-bench,1,280,flap,73.7222\nfolder-r1,synthetic-bench,1,280,gap,4.5038\nfolder-r1,synthetic-bench,1,290,flap,66.0833\nfolder-r1,synthetic-bench,1,290,gap,3.2192\nfolder-r1,synthetic-bench,1,300,flap,58.4444\nfolder-r1,synthetic-bench,1,300,gap,1.9346\nfolder-r1,synthetic-bench,1,310,flap,50.8056\nfolder-r1,synthetic-bench,1,310,gap,0.65\nfolder-r1,synthetic-bench,1,320,flap,43.1667\nfolder-r1,synthetic-bench,1,320,gap,0.65\nfolder-r1,synthetic-bench,1,330,flap,35.5278\nfolder-r1,synthetic-bench,1,330,gap,0.65\nfolder-r1,synthetic-bench,1,340,flap,31.1429\nfolder-r1,synthetic-bench,1,340,gap,0.65\nfolder-r1,synthetic-bench,1,350,flap,27.5714\nfolder-r1,synthetic-bench,1,350,gap,1.825\nfolder-r1,synthetic-bench,1,360,flap,24.0\nfolder-r1,synthetic-bench,1,360,gap,3.0\nfolder-r1,synthetic-bench,1,254,rear_folded,\nfolder-r1,synthetic-bench,1,264,tool_clear,\nfolder-r1,synthetic-bench,2,0,flap,0.0\nfolder-r1,synthetic-bench,2,0,gap,9.0\nfolder-r1,synthetic-bench,2,10,flap,0.0\nfolder-r1,synthetic-bench,2,10,gap,9.0\nfolder-r1,synthetic-bench,2,20,flap,0.0\nfolder-r1,synthetic-bench,2,20,gap,9.0\nfolder-r1,synthetic-bench,2,30,flap,0.0\nfolder-r1,synthetic-bench,2,30,gap,9.0\nfolder-r1,synthetic-bench,2,40,flap,0.0\nfolder-r1,synthetic-bench,2,40,gap,9.0\nfolder-r1,synthetic-bench,2,50,flap,0.0\nfolder-r1,synthetic-bench,2,50,gap,9.0\nfolder-r1,synthetic-bench,2,60,flap,0.0\nfolder-r1,synthetic-bench,2,60,gap,9.0\nfolder-r1,synthetic-bench,2,70,flap,0.0\nfolder-r1,synthetic-bench,2,70,gap,9.0\nfolder-r1,synthetic-bench,2,80,flap,0.0\nfolder-r1,synthetic-bench,2,80,gap,9.0\nfolder-r1,synthetic-bench,2,90,flap,0.0\nfolder-r1,synthetic-bench,2,90,gap,9.0\nfolder-r1,synthetic-bench,2,100,flap,0.0\nfolder-r1,synthetic-bench,2,100,gap,9.0\nfolder-r1,synthetic-bench,2,110,flap,0.0\nfolder-r1,synthetic-bench,2,110,gap,9.0\nfolder-r1,synthetic-bench,2,120,flap,0.0\nfolder-r1,synthetic-bench,2,120,gap,9.0\nfolder-r1,synthetic-bench,2,130,flap,7.4167\nfolder-r1,synthetic-bench,2,130,gap,9.0\nfolder-r1,synthetic-bench,2,140,flap,14.8333\nfolder-r1,synthetic-bench,2,140,gap,9.0\nfolder-r1,synthetic-bench,2,150,flap,22.25\nfolder-r1,synthetic-bench,2,150,gap,9.0\nfolder-r1,synthetic-bench,2,160,flap,29.6667\nfolder-r1,synthetic-bench,2,160,gap,9.0\nfolder-r1,synthetic-bench,2,170,flap,37.0833\nfolder-r1,synthetic-bench,2,170,gap,9.0\nfolder-r1,synthetic-bench,2,180,flap,44.5\nfolder-r1,synthetic-bench,2,180,gap,9.0\nfolder-r1,synthetic-bench,2,190,flap,51.9167\nfolder-r1,synthetic-bench,2,190,gap,9.0\nfolder-r1,synthetic-bench,2,200,flap,59.3333\nfolder-r1,synthetic-bench,2,200,gap,9.0\nfolder-r1,synthetic-bench,2,210,flap,66.75\nfolder-r1,synthetic-bench,2,210,gap,9.0\nfolder-r1,synthetic-bench,2,220,flap,74.1667\nfolder-r1,synthetic-bench,2,220,gap,9.0\nfolder-r1,synthetic-bench,2,230,flap,81.5833\nfolder-r1,synthetic-bench,2,230,gap,9.0\nfolder-r1,synthetic-bench,2,240,flap,89.0\nfolder-r1,synthetic-bench,2,240,gap,9.0\nfolder-r1,synthetic-bench,2,250,flap,89.0\nfolder-r1,synthetic-bench,2,250,gap,8.5385\nfolder-r1,synthetic-bench,2,260,flap,89.0\nfolder-r1,synthetic-bench,2,260,gap,7.6154\nfolder-r1,synthetic-bench,2,270,flap,81.3611\nfolder-r1,synthetic-bench,2,270,gap,6.6923\nfolder-r1,synthetic-bench,2,280,flap,73.7222\nfolder-r1,synthetic-bench,2,280,gap,5.7692\nfolder-r1,synthetic-bench,2,290,flap,66.0833\nfolder-r1,synthetic-bench,2,290,gap,4.8462\nfolder-r1,synthetic-bench,2,300,flap,58.4444\nfolder-r1,synthetic-bench,2,300,gap,3.9231\nfolder-r1,synthetic-bench,2,310,flap,50.8056\nfolder-r1,synthetic-bench,2,310,gap,3.0\nfolder-r1,synthetic-bench,2,320,flap,43.1667\nfolder-r1,synthetic-bench,2,320,gap,3.0\nfolder-r1,synthetic-bench,2,330,flap,35.5278\nfolder-r1,synthetic-bench,2,330,gap,3.0\nfolder-r1,synthetic-bench,2,340,flap,31.1429\nfolder-r1,synthetic-bench,2,340,gap,3.0\nfolder-r1,synthetic-bench,2,350,flap,27.5714\nfolder-r1,synthetic-bench,2,350,gap,3.0\nfolder-r1,synthetic-bench,2,360,flap,24.0\nfolder-r1,synthetic-bench,2,360,gap,3.0\nfolder-r1,synthetic-bench,2,254,rear_folded,\nfolder-r1,synthetic-bench,2,264,tool_clear,\n"
  },
  "reach": {
    "setup": {
      "schema": 1,
      "setup": "reach-r1",
      "name": "Between-frame reach test",
      "source": "synthetic",
      "channels": [
        {
          "id": "gap",
          "label": "Reference clearance",
          "unit": "mm",
          "sense": "min",
          "limit": 0,
          "error": 0.01,
          "rate": 0.01,
          "window": [
            0,
            360
          ]
        }
      ],
      "events": [],
      "handoffs": [],
      "loops": [
        {
          "name": "Narrow-gap regression crank",
          "groundDeg": 0.25,
          "ground": [
            100,
            100
          ],
          "crank": [
            30,
            30
          ],
          "coupler": [
            75,
            75
          ],
          "rocker": [
            54.99994,
            54.99994
          ]
        }
      ]
    },
    "csv": "setup,run,cycle,phase_deg,channel,value\nreach-r1,synthetic-bench,1,0,gap,2\nreach-r1,synthetic-bench,1,10,gap,2\nreach-r1,synthetic-bench,1,20,gap,2\nreach-r1,synthetic-bench,1,30,gap,2\nreach-r1,synthetic-bench,1,40,gap,2\nreach-r1,synthetic-bench,1,50,gap,2\nreach-r1,synthetic-bench,1,60,gap,2\nreach-r1,synthetic-bench,1,70,gap,2\nreach-r1,synthetic-bench,1,80,gap,2\nreach-r1,synthetic-bench,1,90,gap,2\nreach-r1,synthetic-bench,1,100,gap,2\nreach-r1,synthetic-bench,1,110,gap,2\nreach-r1,synthetic-bench,1,120,gap,2\nreach-r1,synthetic-bench,1,130,gap,2\nreach-r1,synthetic-bench,1,140,gap,2\nreach-r1,synthetic-bench,1,150,gap,2\nreach-r1,synthetic-bench,1,160,gap,2\nreach-r1,synthetic-bench,1,170,gap,2\nreach-r1,synthetic-bench,1,180,gap,2\nreach-r1,synthetic-bench,1,190,gap,2\nreach-r1,synthetic-bench,1,200,gap,2\nreach-r1,synthetic-bench,1,210,gap,2\nreach-r1,synthetic-bench,1,220,gap,2\nreach-r1,synthetic-bench,1,230,gap,2\nreach-r1,synthetic-bench,1,240,gap,2\nreach-r1,synthetic-bench,1,250,gap,2\nreach-r1,synthetic-bench,1,260,gap,2\nreach-r1,synthetic-bench,1,270,gap,2\nreach-r1,synthetic-bench,1,280,gap,2\nreach-r1,synthetic-bench,1,290,gap,2\nreach-r1,synthetic-bench,1,300,gap,2\nreach-r1,synthetic-bench,1,310,gap,2\nreach-r1,synthetic-bench,1,320,gap,2\nreach-r1,synthetic-bench,1,330,gap,2\nreach-r1,synthetic-bench,1,340,gap,2\nreach-r1,synthetic-bench,1,350,gap,2\nreach-r1,synthetic-bench,1,360,gap,2\n"
  }
};
