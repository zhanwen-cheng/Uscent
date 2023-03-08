import sys
import getopt
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
from fpdf import FPDF

def usage():
    print('''
          Usage: python3 script.py [option] [patameter]
          -d/--depth          input depth and U counts file
          -o/--outfile        input outfile name
          -h/--help           show possible options
          ''')


opts, args = getopt.getopt(sys.argv[1:], 'hd:o:', ['help', 'depth=', 'outfile='])
for opt, val in opts:
    if opt == '-d' or opt == '--depth':
        depth_file = val
    elif opt == '-o' or opt == '--outfile':
        out_file = val
    elif opt == '-h' or opt == '--help':
        usage()
        sys.exit(1)

####define all need function and expression here
# def df_move1(table,m):
#     table['new_lable']=table.lable
#     table['ori']='#1F78B4'
#     if (max(table.lable)+min(table.lable)) > 2*int(m):
#         table.loc[table.lable<m,'new_lable']= table.loc[table.lable<m,'lable']+max(table.lable)
#         table.loc[table.lable<m,'ori'] = '#DAE3F3'
#     else:
#         table.loc[table.lable>m,'new_lable']= table.loc[table.lable>=m,'lable']-max(table.lable)
#         table.loc[table.lable>m,'ori'] = '#DAE3F3'
#     return(table)
def df_move1(table,m):
    table['new_lable']=table.lable
    table['ori']='#1F78B4'
    table.loc[table.lable>m,'new_lable'] = table.loc[table.lable>m,'lable']-max(table.lable)
    table.loc[table.lable>m,'ori'] = '#DAE3F3'
    return(table)

def df_move2(table,m1,n1,m2,n2):####must sort from min to max
    table['new_lable']=table.lable
    table['ori']='#1F78B4'
    if m2 < n2:
        table.loc[table.lable>m1,'new_lable']= table.loc[table.lable>m1,'lable']-max(table.lable)
        table.loc[table.lable>m1,'ori'] = '#DAE3F3'
    else:
        table.loc[table.lable>n1,'new_lable']= table.loc[table.lable>n1,'lable']-max(table.lable)
        table.loc[table.lable>n1,'ori'] = '#DAE3F3'
    return(table)

def df_move3(table,n,m,k):####must sort from min to max
    table['new_lable']=table.lable
    table['ori']='#1F78B4'
    if (max(table.lable)+min(table.lable)) > 2*int(m):
        table.loc[table.lable<m,'new_lable']= table.loc[table.lable<m,'lable']+max(table.lable)
        table.loc[table.lable<m,'ori'] = '#DAE3F3'
    else:
        table.loc[table.lable>m,'new_lable']= table.loc[table.lable>m,'lable']-max(table.lable)
        table.loc[table.lable>m,'ori'] = '#DAE3F3'
    return(table)

#def func1(x,h,a,u,sig):
#    return h+(a*np.exp(-(x - u) ** 2 / (2 * sig ** 2)))/ (sig * math.sqrt(2 * math.pi))
def func2(x,u0,sig0,a,u,sig):
    return (np.exp(-(x - u0) ** 2 / (2 * sig0 ** 2)))/ (sig0 * math.sqrt(2 * math.pi))+(a*np.exp(-(x - u) ** 2 / (2 * sig ** 2)))/(sig * math.sqrt(2 * math.pi))


################################################################################
###############################first page for poly fit##########################
################################################################################
pdf=FPDF(format='letter', unit='in')
pdf.add_page()
df_new = pd.read_table(depth_file, sep=',')
poly_4 = np.polyfit(df_new.lable,df_new.Bcounts,4)
poly_4_item = np.poly1d(poly_4)

####insert table in PDF####
# Effective page width, or just epw
epw = pdf.w - 2*pdf.l_margin
# Set column width to 1/4 of effective page width to distribute content
# evenly across table and page
col_width = epw/5
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The parameter of quartic polynomial', align='C')
pdf.ln(0.2)
pdf.image('/data/home/zhanwen/software/script/fx1.png',w=8,h=0.8)
pdf.ln(0.2)
pdf.set_font('Times','',8.0)
data1 = [['The quartic item X4','The cubic item X3','The quadratic iterm X2','The primary item X','The constant'],
[poly_4_item[4],poly_4_item[3],poly_4_item[2],poly_4_item[1],poly_4_item[0]]]
#print(data1)
th = pdf.font_size
for row in data1:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)
pdf.ln(0.2)

yvals=poly_4_item(df_new.lable)
max_y = max(yvals)
min_y = min(yvals)

der_poly_4_item=list(np.polyder(poly_4_item,1).r)
der_poly_4_item=[x for x in der_poly_4_item if x >=0]
der_poly_4_item=[x for x in der_poly_4_item if x <= max(df_new.lable)]
der_poly_4_item.sort()

if isinstance(der_poly_4_item[0],complex):
    Xlist=list(np.round([x.real for x in der_poly_4_item if (x.imag==0)],2))
    Ylist=list(poly_4_item(Xlist))
else:
    Xlist=list(np.round([x.real for x in der_poly_4_item],2))
    Ylist=list(poly_4_item(Xlist))

if len(Xlist)==1:
    df_new=df_move1(df_new,Xlist[0])
elif len(Xlist)==2:
    df_new=df_move2(df_new,Xlist[0],Xlist[1],Ylist[0],Ylist[1])
elif len(Xlist)==3:
    df_new=df_move3(df_new,Xlist[0],Xlist[1],Xlist[2])

plt.figure(figsize=(6,4),dpi=1000)
plt.scatter(df_new.lable,df_new.Bcounts,c=df_new.ori,s=4)
plt.scatter(Xlist,Ylist,marker='^',color='r')
plt.plot(df_new.lable,yvals,'r',label='poly fit')
plt.savefig('./Bcounts_poly4.jpg')
pdf.image('./Bcounts_poly4.jpg',w=6, h=4)
pdf.ln(0.2)

data2=[
        ['X'+str(i) for i in range(1,len(Xlist)+1)],Xlist,
        ['Y'+str(i) for i in range(1,len(Ylist)+1)],Ylist
        ]

pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'Where the derivative is zero', align='C')
pdf.ln(0.2)
pdf.set_font('Times','',8.0)
th = pdf.font_size
for row in data2:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)
pdf.ln(0.2)


################################################################################
#############################another page for model fit#########################
################################################################################
df_new=df_new.sort_values(by='new_lable',ascending=True).reset_index(drop=True)
#popthunt1, pcovhunt1 = curve_fit(func1, df_new.new_lable, df_new.Bcounts,maxfev=50000) 
popthunt2, pcovhunt2 = curve_fit(func2, df_new.new_lable, df_new.Bcounts,maxfev=50000) 

#hhunt1 = popthunt1[0]
#ahunt1 = popthunt1[1]
#uhunt1 = popthunt1[2]
#sighunt1 = popthunt1[3]


u0hunt2 = popthunt2[0]
sig0hunt2 = popthunt2[1]
ahunt2 = popthunt2[2]
uhunt2 = popthunt2[3]
sighunt2 = popthunt2[4]


#yhuntvals1 = func1(df_new.new_lable,hhunt1,ahunt1,uhunt1,sighunt1)
yhuntvals2 = func2(df_new.new_lable,u0hunt2,sig0hunt2,ahunt2,uhunt2,sighunt2)

index1 = (ahunt2*np.exp(-1/(2 * sighunt2 ** 2)) / (sighunt2 * math.sqrt(2 * math.pi)))/(max_y-min_y)
index2 = (ahunt2*1/sighunt2)/(max_y-min_y)
index3 = (ahunt2*2/(np.exp(sighunt2)-1))/(max_y-min_y)
index4 = (ahunt2*1/(pow(sighunt2,1/3)))/(max_y-min_y)
index5 = (ahunt2*2/(np.exp(math.sqrt(sighunt2))-1))/(max_y-min_y)

pdf.add_page()
data3=[
        ['u0','sig0','a1','u1','sig1'],
        [popthunt2[0],popthunt2[1],popthunt2[2],popthunt2[3],popthunt2[4]],
        ]
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The parameter of normal distribution', align='C')
pdf.ln(0.2)
pdf.image('/data/home/zhanwen/software/script/fx2.png',w=8,h=0.8)
pdf.ln(0.2)
epw = pdf.w - 2*pdf.l_margin
col_width = epw/5
pdf.set_font('Times','',8.0)
th = pdf.font_size
for row in data3:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)
pdf.ln(0.2)


plt.figure(figsize=(6,4),dpi=1000)
plt.scatter(df_new.new_lable,df_new.Bcounts,c=df_new.ori,s=4)
plt.plot(df_new.new_lable,yhuntvals2,'r',label='normal fit')
plt.savefig('./Bcounts_fit_pic.jpg')
pdf.image('./Bcounts_fit_pic.jpg',w=6, h=4)
pdf.ln(0.2)

data4=[
        ['index1','index2','index3','index4','index5'],
        [index1,index2,index3,index4,index5],
        ]
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The index calculated from normal distribution', align='C')
pdf.ln(0.2)
pdf.image('/data/home/zhanwen/software/script/index.png',w=8,h=1.5)
pdf.ln(0.2)
epw = pdf.w - 2*pdf.l_margin
col_width = epw/5
pdf.set_font('Times','',8.0)
th = pdf.font_size
for row in data4:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)


pdf.output(out_file)

































