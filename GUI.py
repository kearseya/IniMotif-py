if int(numberofrunsinput.get()) == len(namesindirectory):
                filenamesinputs.delete(0,10)
                filenamesinputs.insert(0, firstfileauto)
                for x in range(0, (int(numberofrunsinput.get())-1)):
                    for j in namesindirectory:
                        threeletter = str(firstfileauto[:3])
                        sixnumbers = str(int(firstfileauto[3:-6])-(x+1))
                        if len(sixnumbers) < 6:
                            sixnumbers = (6-len(sixnumbers))*"0"+sixnumbers
                        namenoext = threeletter+sixnumbers
                        if namenoext in j[:-6]:
                            work = SELEXlist[x*2]
                            work.delete(0,10)
                            work.insert(0, j)
            if int(numberofrunsinput.get()) < len(namesindirectory):
                filenamesinputs.delete(0,10)
                threeletterstart = str(firstfileauto[:3])
                sixnumbersstart = str(int(firstfileauto[3:-6])-(int(startroundinputs.get())-1))
                if len(sixnumbersstart) < 6:
                    sixnumbersstart = (6-len(sixnumbersstart))*"0"+sixnumbersstart
                startingfile = threeletterstart+sixnumbersstart+firstfileauto[-6:]
                filenamesinputs.insert(0, startingfile)
                for x in range(0, int(numberofrunsinput.get())-1):
                    for j in namesindirectory:
                        threeletter = str(startingfile[:3])
                        sixnumbers = str(int(startingfile[3:-6])-(x+1))
                        if len(sixnumbers) < 6:
                            sixnumbers = (6-len(sixnumbers))*"0"+sixnumbers
                        namenoext = threeletter+sixnumbers
                        if namenoext in j[:-6]:
                            work = SELEXlist[x*2]
                            work.delete(0,10)
                            work.insert(0, j)
