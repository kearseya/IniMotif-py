from tkinter import *


window = Tk()

titleframe = Frame(window)
titleframe.pack(side=TOP)

title = Label(titleframe, text="IniMotif-py")
title.pack(side=TOP)

usage = 0

initialdetailsframe = Frame(window)
initialdetailsframe.pack(side=TOP, anchor="w")

variableinputform = Frame(window)
variableinputform.pack()


experimenttype = Label(initialdetailsframe, text="Experiment type: ")
experimenttype.grid(row=7, column=0, sticky="e")
experimentformatlist = ["SELEX", "ChIP", "ATAC", "DNase"]
valueforexdrop = StringVar()
valueforexdrop.set(experimentformatlist[0])
experimenttypedrop = OptionMenu(initialdetailsframe, valueforexdrop, "SELEX", "ChIP", "ATAC", "DNase")
experimenttypedrop.grid(row=7, column=1)


filenameslabels = Label(variableinputform, text="File name(s): ", anchor="w")
filenameslabels.grid(row=0, column=1)
filenamesinputs = Entry(variableinputform, textvariable=StringVar())
filenamesinputs.grid(row=1, column=1)

readlengthslabels = Label(variableinputform, text="Read lengths: ")
readlengthslabels.grid(row=0, column=2)
readlengthsinputs = Entry(variableinputform, textvariable=IntVar())
readlengthsinputs.delete(0)
readlengthsinputs.grid(row=1, column=2)

fiveprimebarlabels = Label(variableinputform, text="5' barcode/primer: ")
fiveprimebarlabels.grid(row=0, column=3)
fiveprimebarinputs = Entry(variableinputform, textvariable=StringVar())
fiveprimebarinputs.delete(0)
fiveprimebarinputs.grid(row=1, column=3)

threeprimebarlabels = Label(variableinputform, text="3' barcode/primer: ")
threeprimebarlabels.grid(row=0, column=4)
threeprimebarinputs = Entry(variableinputform, textvariable=StringVar())
threeprimebarinputs.delete(0)
threeprimebarinputs.grid(row=1, column=4)

multiroundsamefilelabels = Label(initialdetailsframe, text="Multiple rounds in same file: ")
multiroundsamefilelabels.grid(row=8, column=0, sticky="e", padx=(20,0))
multiround = BooleanVar()
multiroundcheckbox = Checkbutton(initialdetailsframe, variable=multiround)
multiroundcheckbox.grid(row=8, column=1)
multiroundcheckbox.select()

knownbarcodeslabels = Label(initialdetailsframe, text="Known barcode: ")
knownbarcodeslabels.grid(row=9, column=0, sticky="e", padx=(20,0))
knownbarcodes = BooleanVar()
knownbarcodescheckbox = Checkbutton(initialdetailsframe, variable=knownbarcodes)
knownbarcodescheckbox.grid(row=9, column=1)
knownbarcodescheckbox.select()


def add_rows():
    global SELEXlist
    global fileinputslist
    global fiveprimebar
    global threeprimebar
    x = int(numberofrunsinput.get())
    SELEXlist = []
    fileinputslist = []
    fiveprimebar = []
    threeprimebar = []
    extype = valueforexdrop.get()
    knownbarcodesval = knownbarcodes.get()
    #global filenameslabels
    #global filenamesinputs
    #if extype == "SELEX":
        #global readlengthslabels
        #global readlengthsinputs
    if extype == "SELEX":
        for x in range(0, (x-1)*2):
            SELEXlist.append("position"+str(x))
    else:
        for x in range(0,x-1):
            fileinputslist.append("file"+str(x))
    x = int(numberofrunsinput.get())
    if knownbarcodesval == True:
        for x in range(0, x-1):
            fiveprimebar.append("5barcode"+str(x))
            threeprimebar.append("3barcode"+str(x))
    x = int(numberofrunsinput.get())
    if extype == "SELEX":
        for x in range(0, (x-1)*2):
            if x%2 == 0:
                SELEXlist[x] = Entry(variableinputform, textvariable=StringVar())
                if usage%2 == 0:
                    SELEXlist[x].delete(0, 10)
                SELEXlist[x].grid(row=(x//2)+2, column=(x%2)+1)
            if x%2 == 1:
                SELEXlist[x] = Entry(variableinputform, textvariable=IntVar())
                if usage%2 == 0:
                    SELEXlist[x].delete(0, 10)
                SELEXlist[x].grid(row=(x//2)+2, column=(x%2)+1)
    else:
        for x in range(0, x-1):
            fileinputslist[x] = Entry(variableinputform, textvariable=StringVar())
            if usage%2 == 0:
                fileinputslist[x].delete(0, 10)
            fileinputslist[x].grid(row=x+2, column=1)
            readlengthslabels.destroy()
            readlengthsinputs.destroy()
    if knownbarcodesval == True:
        for x in range(0, x-1):
            fiveprimebar[x] = Entry(variableinputform, textvariable=StringVar())
            threeprimebar[x] = Entry(variableinputform, textvariable=StringVar())
            if usage%2 == 0:
                fiveprimebar[x].delete(0, 20)
                threeprimebar[x].delete(0, 20)
            fiveprimebar[x].grid(row=x+2, column=3)
            threeprimebar[x].grid(row=x+2, column=4)
    if knownbarcodesval == False:
        fiveprimebarlabels.destroy()
        fiveprimebarinputs.destroy()
        threeprimebarlabels.destroy()
        threeprimebarinputs.destroy()



identifiernamelabel = Label(initialdetailsframe, text="Analysis identifier: ")
identifiernamelabel.grid(row=0, column=0, sticky="e")
identifiernameinput = Entry(initialdetailsframe)
identifiernameinput.grid(row=0, column=1)

pathtodirectorylabel = Label(initialdetailsframe, text="Path to directory: ")
pathtodirectorylabel.grid(row=1, column=0, sticky="e")
pathtodirectoryinput = Entry(initialdetailsframe)
pathtodirectoryinput.grid(row=1, column=1)

numberofrunslabel = Label(initialdetailsframe, text="Number of runs: ")
numberofrunslabel.grid(row=2, column=0, sticky="e")
numberofrunsinput = Entry(initialdetailsframe, textvariable=IntVar())
numberofrunsinput.delete(0)
numberofrunsinput.grid(row=2, column=1)

reversecomplementwantedlabel = Label(initialdetailsframe, text="Reverse compliment: ")
reversecomplementwantedlabel.grid(row=3, column=0, sticky="e", padx=(20,0))
reversecomplementwanted = BooleanVar()
reversecomplementwantedcheckbox = Checkbutton(initialdetailsframe, variable=reversecomplementwanted)
reversecomplementwantedcheckbox.grid(row=3, column=1)

logoformatlabel = Label(initialdetailsframe, text="Logo format: ")
logoformatlabel.grid(row=4, column=0, sticky="e")
logoformatlist = ["bits", "frequency"]
valueforlogodrop = StringVar()
valueforlogodrop.set(logoformatlist[0])
logoformatdropdown = OptionMenu(initialdetailsframe, valueforlogodrop, "bits", "frequency")
logoformatdropdown.grid(row=4, column=1)

formchangebutton = Button(initialdetailsframe, text="Enter", pady=1, command=add_rows)
formchangebutton.grid(row=2, column=2, padx=10)

inimotifimage = PhotoImage(file='figures/GUIgraphics/tobylogo.png')
inimotifimagelabel = Label(initialdetailsframe, image=inimotifimage)
inimotifimagelabel.grid(row=1, column=3, padx=20)

minimumkvalueslabels = Label(initialdetailsframe, text="Min K: ")
minimumkvalueslabels.grid(row=5, column=0,  sticky="e")
minimimkvaluesinputs = Entry(initialdetailsframe, textvariable=IntVar())
minimimkvaluesinputs.delete(0)
minimimkvaluesinputs.grid(row=5, column=1)

maximumkvalueslabels = Label(initialdetailsframe, text="Max K: ")
maximumkvalueslabels.grid(row=6, column=0,  sticky="e")
maximumkvaluesinputs = Entry(initialdetailsframe, textvariable=IntVar())
maximumkvaluesinputs.delete(0)
maximumkvaluesinputs.grid(row=6, column=1)

startroundlabels = Label(initialdetailsframe, text="Start round: ")
startroundlabels.grid(row=10, column=0,  sticky="e")
startroundinputs = Entry(initialdetailsframe, textvariable=IntVar())
startroundinputs.delete(0)
startroundinputs.insert(0, int(1))
startroundinputs.grid(row=10, column=1)


"""
filenameslabels = Label(variableinputform, text="File name: ", anchor="w")
filenameslabels.grid(row=0, column=1)
filenamesinputs = Entry(variableinputform, textvariable=StringVar())
filenamesinputs.grid(row=1, column=1)
"""
"""
readlengthslabels = Label(variableinputform, text="Read lengths: ")
readlengthslabels.grid(row=1, column=2)
readlengthsinputs = Entry(variableinputform, textvariable=IntVar())
readlengthsinputs.delete(0)
readlengthsinputs.grid(row=2, column=2)
"""



def autofiller():
    usage =+ 1
    extype = valueforexdrop.get()
    #if usage%2 != 0:
        #add_rows()
    if int(numberofrunsinput.get()) != 0:
        firstfile = str(filenamesinputs.get())
        if extype == "SELEX":
            readval = int(readlengthsinputs.get())
            for i in range(0, (int(numberofrunsinput.get())-1)*2):
                if i%2 == 0:
                    work = SELEXlist[i]
                    if firstfile[-1].isdigit() == True:
                        endnum = int(firstfile[-1])+(i//2)+1
                        predicted = str(firstfile[:-1])+str(endnum)
                    work.insert(0, string=predicted)
                if i%2 == 1:
                    work = SELEXlist[i]
                    work.insert(0, string=readval)
        else:
            for i in range(0, (int(numberofrunsinput.get())-1)):
                work = fileinputslist[i]
                if firstfile[-1].isdigit() == True:
                    endnum = int(firstfile[-1])+i+1
                    predicted = str(firstfile[:-1])+str(endnum)
                    work.insert(0, string=predicted)


def makeorderedinputlist():
    global inputlist
    inputlist = []

    inputlist.append(str(identifiernameinput.get()))
    inputlist.append(str(pathtodirectoryinput.get()))
    inputlist.append(int(numberofrunsinput.get()))
    inputlist.append(bool(reversecomplementwanted.get()))
    #inputlist.append(str(filenamesinputs.get()))
    #inputlist.append(int(runnumbersinputs.get()))
    inputlist.append(int(minimimkvaluesinputs.get()))
    inputlist.append(int(maximumkvaluesinputs.get()))
    inputlist.append(str(valueforexdrop.get()))
    inputlist.append(str(filenamesinputs.get()))
    extype = str(valueforexdrop.get())
    if extype == "SELEX":
        inputlist.append(int(readlengthsinputs.get()))
        for x in range(0, (int(numberofrunsinput.get())-1)*2):
            position = SELEXlist[x]
            if x%2 == 0:
                inputlist.append(str(position.get()))
            if x%2 == 1:
                try:
                    inputlist.append(int(position.get()))
                except:
                    continue
    else:
        for x in range(0, int(numberofrunsinput.get())-1):
            position = fileinputslist[x]
            inputlist.append(str(position.get()))

    global logotype
    logotype = str(valueforlogodrop.get())

    global startround
    startround = int(startroundinputs.get())
    global multiround
    multiround = multiround.get()
    global knownbarcode
    knownbarcode = knownbarcodes.get()

    if knownbarcode == True:
        inputlist.append(str(fiveprimebarinputs.get()))
        inputlist.append(str(threeprimebarinputs.get()))
        for x in range(0, x-1):
            inputlist.append(fiveprimebar[x].get())
            inputlist.append(threeprimebar[x].get())

    return inputlist


#autofilldetailslabel = Label(window, text="Autofill: ")
#autofilldetailslabel.pack(side=LEFT, anchor="n", text="Autofill")
autofilldetailscheckbox = Button(window, command=autofiller, text="Autofill")
autofilldetailscheckbox.pack(side=LEFT, anchor="n", padx=20, pady=20)


submitdetailscheckbox = Button(window, command=makeorderedinputlist, text="SUBMIT")
submitdetailscheckbox.pack(side=RIGHT, anchor="n", padx=20, pady=20)


window.mainloop()




"""
def on_entry_click(event):
    #function that gets called whenever entry is clicked#
    if entry.get() == 'Enter your user name...':
       entry.delete(0, "end") # delete all the text in the entry
       entry.insert(0, '') #Insert blank for user input
       entry.config(fg = 'black')
def on_focusout(event):
    if entry.get() == '':
        entry.insert(0, 'Enter your username...')
        entry.config(fg = 'grey')
"""
