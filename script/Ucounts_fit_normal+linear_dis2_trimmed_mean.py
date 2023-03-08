import sys
import getopt
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
from fpdf import FPDF
from scipy import stats

def usage():
    print('''
          Usage: python3 script.py [option] [patameter]
          -d/--depth          input depth and U counts file
          -g/--gc             input gc file
          -o/--outfile        input outfile name
          -h/--help           show possible options
          ''')


opts, args = getopt.getopt(sys.argv[1:], 'hd:g:o:', ['help', 'depth=', 'gc=','outfile='])
for opt, val in opts:
    if opt == '-d' or opt == '--depth':
        depth_file = val
    elif opt == '-o' or opt == '--outfile':
        out_file = val
    elif opt == '-g' or opt == '--gc':
        gc_file = val
    elif opt == '-h' or opt == '--help':
        usage()
        sys.exit(1)

def df_move1(table,m):
    table['new_lable']=table.lable
    table['ori']='#1F78B4'
    table.loc[table.lable>m,'new_lable'] = table.loc[table.lable>m,'lable']-max(table.lable)
    table.loc[table.lable>m,'ori'] = '#DAE3F3'
    der_p=max(table.lable)-max(table.lable)
    return(table,der_p)

def df_move2(table,m1,n1,m2,n2):####must sort from min to max
    table['new_lable']=table.lable
    table['ori']='#1F78B4'
    if m2 < n2:
        table.loc[table.lable>m1,'new_lable']= table.loc[table.lable>m1,'lable']-max(table.lable)
        table.loc[table.lable>m1,'ori'] = '#DAE3F3'
        der_p=n1-max(table.lable)
    else:
        table.loc[table.lable>n1,'new_lable']= table.loc[table.lable>n1,'lable']-max(table.lable)
        table.loc[table.lable>n1,'ori'] = '#DAE3F3'
        der_p=m1
    return(table,der_p)

def df_move3(table,n1,m1,k1,n2,m2,k2):####must sort from min to max
    table['new_lable']=table.lable
    table['ori']='#1F78B4'
    if max(n2,m2,k2)==n2:
        table.loc[table.lable>m1,'new_lable']= table.loc[table.lable>m1,'lable']-max(table.lable)
        table.loc[table.lable>m1,'ori'] = '#DAE3F3'
        der_p=n1
    elif max(n2,m2,k2)==k2:
        table.loc[table.lable<m1,'new_lable']= table.loc[table.lable<m1,'lable']+max(table.lable)
        table.loc[table.lable<m1,'ori'] = '#DAE3F3'
        der_p=k1
    else:
        if n2>k2:
            table.loc[table.lable>k1,'new_lable']= table.loc[table.lable>k1,'lable']-max(table.lable)
            table.loc[table.lable>k1,'ori'] = '#DAE3F3'
            der_p=m1
        elif n2<k2:
            table.loc[table.lable<n1,'new_lable']= table.loc[table.lable<n1,'lable']+max(table.lable)
            table.loc[table.lable<n1,'ori'] = '#DAE3F3'
            der_p=m1
    return(table,der_p)

def func1(x,a,b):
    return (a*x+b)
def func2(x,a,u,sig):
    return (a*np.exp(-(x - u) ** 2 / (2 * sig ** 2)))/(sig * math.sqrt(2 * math.pi))


################################################################################
###############################first page for poly fit##########################
################################################################################
pdf=FPDF(format='letter', unit='in')
pdf.add_page()
df_new = pd.read_table(depth_file, sep=',')
df_new = df_new[0:-1]
gc_table = pd.read_table(gc_file,sep=",")
gc_table.gc = gc_table.gc/100
df_new['ori_Ucounts']=df_new.Ucounts
df_new['Ucounts'] = df_new.Ucounts*gc_table.gc
df_new.Ucounts=df_new.Ucounts/(stats.trim_mean(df_new.ori_Ucounts,0.05))
df_new['ori_lable']=df_new.lable
df_new.lable=df_new.lable/(stats.trim_mean(df_new.ori_lable,0.05))
poly_4 = np.polyfit(df_new.lable,df_new.Ucounts,4)
poly_4_item = np.poly1d(poly_4)

####insert table in PDF####
# Effective page width, or just epw
epw = pdf.w - 2*pdf.l_margin
# Set column width to 1/4 of effective page width to distribute content
# evenly across table and page
col_width = epw/5

# insert the original dirtribution plot
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The original distribution scatter plot', align='C')
pdf.ln(0.2)
plt.figure(figsize=(5,2),dpi=800)
plt.scatter(df_new.ori_lable,df_new.ori_Ucounts,s=2)
plt.savefig('./Ucounts.jpg')
pdf.image('./Ucounts.jpg',w=5, h=2)
pdf.ln(0.2)

# insert polyfit
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
    df_new_data=df_move1(df_new,Xlist[0])
    df_new=df_new_data[0]
    der_point=df_new_data[1]
elif len(Xlist)==2:
    df_new_data=df_move2(df_new,Xlist[0],Xlist[1],Ylist[0],Ylist[1])
    df_new=df_new_data[0]
    der_point=df_new_data[1]
elif len(Xlist)==3:
    df_new_data=df_move3(df_new,Xlist[0],Xlist[1],Xlist[2],Ylist[0],Ylist[1],Ylist[2])
    df_new=df_new_data[0]
    der_point=df_new_data[1]

plt.figure(figsize=(6,4),dpi=1000)
plt.scatter(df_new.lable,df_new.Ucounts,c=df_new.ori,s=4)
plt.scatter(Xlist,Ylist,marker='^',color='r')
plt.plot(df_new.lable,yvals,'r',label='poly fit')
plt.savefig('./Ucounts_poly4.jpg')
pdf.image('./Ucounts_poly4.jpg',w=6, h=4)
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
##########################second page for linear model fit######################
################################################################################
df_new=df_new.sort_values(by='new_lable',ascending=True).reset_index(drop=True)
df_new_left=df_new.loc[df_new.new_lable<=der_point,]
df_new_right=df_new.loc[df_new.new_lable>=der_point,]
popthunt_left, pcovhunt_left = curve_fit(func1, df_new_left.new_lable, df_new_left.Ucounts,maxfev=50000)
popthunt_right, pcovhunt_right = curve_fit(func1, df_new_right.new_lable, df_new_right.Ucounts,maxfev=50000)


a_left=popthunt_left[0]
b_left=popthunt_left[1]
a_right=popthunt_right[0]
b_right=popthunt_right[1]

yhuntvals_left = func1(df_new_left.new_lable,a_left,b_left)
yhuntvals_right = func1(df_new_right.new_lable,a_right,b_right)


pdf.add_page()
data3=[
        ['a_left','b_left','a_right','b_right'],
        [a_left,b_left,a_right,b_right],
        ]
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The parameter of linear distribution', align='C')
pdf.ln(0.2)
pdf.image('/data/home/zhanwen/software/script/fx3.png',w=8,h=0.8)
pdf.ln(0.2)
epw = pdf.w - 2*pdf.l_margin
col_width = epw/4
pdf.set_font('Times','',8.0)
th = pdf.font_size
for row in data3:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)
pdf.ln(0.2)


plt.figure(figsize=(6,4),dpi=1000)
plt.scatter(df_new.new_lable,df_new.Ucounts,c=df_new.ori,s=4)
plt.plot(df_new_left.new_lable,yhuntvals_left,'r',label='left linear fit')
plt.plot(df_new_right.new_lable,yhuntvals_right,'r',label='right linear fit')
plt.savefig('./Ucounts_linear_fit_pic.jpg')
pdf.image('./Ucounts_linear_fit_pic.jpg',w=6, h=4)
pdf.ln(0.2)

data4=[
        ['index1','index2','index3'],
        [abs(a_left),abs(a_right),(abs(a_left)+abs(a_right))/2],
        ]
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The index calculated from linear distribution', align='C')
pdf.ln(0.2)
#pdf.image('/data/home/zhanwen/software/script/index.png',w=8,h=1.5)
#pdf.ln(0.2)
epw = pdf.w - 2*pdf.l_margin
col_width = epw/3
pdf.set_font('Times','',8.0)
th = pdf.font_size
for row in data4:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)

################################################################################
##########################third page for noraml model fit#######################
################################################################################
param_bounds=([0,-np.inf,0],
              [np.inf,np.inf,np.inf])
popthunt2, pcovhunt2 = curve_fit(func2, df_new.new_lable, df_new.Ucounts,bounds=param_bounds,maxfev=50000) 

ahunt = popthunt2[0]
uhunt = popthunt2[1]
sighunt = popthunt2[2]
#print("a:",ahunt)
#print("u:",uhunt)
#print("sig:",sighunt)

yhuntvals2 = func2(df_new.new_lable,ahunt,uhunt,sighunt)

index1 = (ahunt*1/sighunt)
index2 = (ahunt*2/(np.exp(sighunt)-1))
index3 = (ahunt*1/(pow(sighunt,1/3)))
index4 = (ahunt*2/(np.exp(math.sqrt(sighunt))-1))

pdf.add_page()
data5=[
        ['a','u','sig'],
        [popthunt2[0],popthunt2[1],popthunt2[2]],
        ]
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The parameter of normal distribution', align='C')
pdf.ln(0.2)
pdf.image('/data/home/zhanwen/software/script/fx2.png',w=8,h=0.8)
pdf.ln(0.2)
epw = pdf.w - 2*pdf.l_margin
col_width = epw/3
pdf.set_font('Times','',8.0)
th = pdf.font_size
for row in data5:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)
pdf.ln(0.2)


plt.figure(figsize=(6,4),dpi=1000)
plt.scatter(df_new.new_lable,df_new.Ucounts,c=df_new.ori,s=4)
plt.plot(df_new.new_lable,yhuntvals2,'r',label='normal fit')
plt.savefig('./Ucounts_normal_fit_pic.jpg')
pdf.image('./Ucounts_normal_fit_pic.jpg',w=6, h=4)
pdf.ln(0.2)

data6=[
        ['index1','index2','index3','index4'],
        [index1,index2,index3,index4],
        ]
pdf.set_font('Times','B',14.0)
pdf.cell(epw, 0.0, 'The index calculated from normal distribution', align='C')
pdf.ln(0.2)
pdf.image('/data/home/zhanwen/software/script/index.png',w=8,h=1.5)
pdf.ln(0.2)
epw = pdf.w - 2*pdf.l_margin
col_width = epw/4
pdf.set_font('Times','',8.0)
th = pdf.font_size
for row in data6:
    for datum in row:
        pdf.cell(col_width, 2*th, str(datum), border=1, align='C')
    pdf.ln(2*th)


pdf.output(out_file)

































