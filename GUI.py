from tkinter import *


usage = 0

def add_rows():
    global efflist
    x = int(numberofrunsinput.get())
    efflist = []
    for x in range(0,x+1):
        for j in range(0,5):
            efflist.append("entry"+str(x+2)+"fields"+str(j))
    for x in range(0,x-1):
        for j in range(0,5):
            if j == 0:
                efflist[(x*5)+j] = Entry(variableinputform, textvariable=StringVar())
                if usage%2 == 0:
                    efflist[(x*5)+j].delete(0, 10)
                efflist[(x*5)+j].grid(row=x+2, column=j)
            if j > 0:
                efflist[(x*5)+j] = Entry(variableinputform, textvariable=IntVar())
                if usage%2 == 0:
                    efflist[(x*5)+j].delete(0, 10)
                efflist[(x*5)+j].grid(row=x+2, column=j)
            if j == 1:
                efflist[(x*5)+j].insert(string=str(x+2), index=1)






window = Tk()

titleframe = Frame(window)
titleframe.pack(side=TOP)

title = Label(titleframe, text="IniMotif-py")
title.pack(side=TOP)


initialdetailsframe = Frame(window)
initialdetailsframe.pack(side=TOP, anchor="w")


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

formchangebutton = Button(initialdetailsframe, text="Enter", pady=1, command=add_rows)
formchangebutton.grid(row=2, column=2, padx=10)

inimotifimage = PhotoImage(file='figures/GUIgraphics/tobylogo.png')
inimotifimagelabel = Label(initialdetailsframe, image=inimotifimage)
inimotifimagelabel.grid(row=1, column=3, padx=20)


variableinputform = Frame(window)
variableinputform.pack()


filenameslabels = Label(variableinputform, text="File name: ", anchor="w")
filenameslabels.grid(row=0, column=0)
filenamesinputs = Entry(variableinputform, textvariable=StringVar())
filenamesinputs.grid(row=1, column=0)

runnumberslabels = Label(variableinputform, text="Run number: ")
runnumberslabels.grid(row=0, column=1)
runnumbersinputs = Entry(variableinputform, textvariable=IntVar())
runnumbersinputs.delete(0)
runnumbersinputs.insert(0, 1)
runnumbersinputs.grid(row=1, column=1)

readlengthslabels = Label(variableinputform, text="Read lengths: ")
readlengthslabels.grid(row=0, column=2)
readlengthsinputs = Entry(variableinputform, textvariable=IntVar())
readlengthsinputs.delete(0)
readlengthsinputs.grid(row=1, column=2)

minimumkvalueslabels = Label(variableinputform, text="Min K: ")
minimumkvalueslabels.grid(row=0, column=3)
minimimkvaluesinputs = Entry(variableinputform, textvariable=IntVar())
minimimkvaluesinputs.delete(0)
minimimkvaluesinputs.grid(row=1, column=3)

maximumkvalueslabels = Label(variableinputform, text="Max K: ")
maximumkvalueslabels.grid(row=0, column=4)
maximumkvaluesinputs = Entry(variableinputform, textvariable=IntVar())
maximumkvaluesinputs.delete(0)
maximumkvaluesinputs.grid(row=1, column=4)


def autofiller():
    usage =+ 1
    if usage%2 != 0:
        add_rows()
    if int(numberofrunsinput.get()) != 0:
        readval = int(readlengthsinputs.get())
        minkvalue = int(minimimkvaluesinputs.get())
        maxkvalue = int(maximumkvaluesinputs.get())
        for i in range(0, (int(numberofrunsinput.get())+1)):
            for j in range(0,5):
                if j == 2:
                    work = efflist[(i*5)+j]
                    work.insert(0, string=readval)
                if j == 3:
                    work = efflist[(i*5)+j]
                    work.insert(0, string=minkvalue)
                if j == 4:
                    work = efflist[(i*5)+j]
                    work.insert(0, string=maxkvalue)


def makeorderedinputlist():
    global inputlist
    inputlist = []
    inputlist.append(str(identifiernameinput.get()))
    inputlist.append(str(pathtodirectoryinput.get()))
    inputlist.append(int(numberofrunsinput.get()))
    inputlist.append(str(reversecomplementwanted.get()))
    inputlist.append(str(filenamesinputs.get()))
    inputlist.append(int(runnumbersinputs.get()))
    inputlist.append(int(readlengthsinputs.get()))
    inputlist.append(int(minimimkvaluesinputs.get()))
    inputlist.append(int(maximumkvaluesinputs.get()))
    for x in range(0, int(numberofrunsinput.get())-1):
        for j in range(0,5):
            position = efflist[(x*5)+j]
            if j == 0:
                inputlist.append(str(position.get()))
            if j > 0:
                inputlist.append(int(position.get()))
    print(inputlist)
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