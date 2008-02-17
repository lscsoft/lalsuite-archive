#!/usr/bin/env python

import copy
import sys
import os.path
import string
import math
from PyQt4 import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import powerflux
from powerflux import compute_scores
from powerflux import compute_matched_snr
from powerflux import compute_single_snr
from powerflux import get_power_sum

class candidate:
	def __init__(self):
		pass
	def dump(self):
		def item_print(item):
			a=getattr(self, item)
			print item + "\t"+str(a)
		map(item_print, dir(self));
		
	def coords(self):
		return (self.frequency, self.spindown, self.ra, self.dec, self.iota, self.psi)
	
	def copy_coords(self, c):
		self.frequency=c.frequency
		self.spindown=c.spindown
		self.ra=c.ra
		self.dec=c.dec
		self.iota=c.iota
		self.psi=c.psi

	
value_cache={}

value_mode="matched"

def compute(c, mode="default"):
	if mode=="default" : mode=value_mode
	if mode=="full" :
		mode=value_mode
		if mode=="matched" :
			mode="full"
	
	ans=value_cache.get((mode, c.coords()))
	if(ans==None):
		ans=copy.copy(c)
		if mode=="single":
			compute_single_snr(ans)
		elif mode=="matched":
			compute_matched_snr(ans)
		elif mode=="full":
			compute_scores(ans)
		else :
			print "Unknown mode"+value_mode
		value_cache[(mode, c.coords())]=ans
	return ans

def set_value_mode(mode):
	global value_mode
	value_mode=mode

def set_veto(start_gps, stop_gps, flag, dataset="ALL"):
	if dataset=="ALL":
		dataset= -1
	else:
		for i in range(len(datasets)):
			if datasets[i]["name"]==dataset:
				dataset=i
				break
	a=powerflux.set_veto(start_gps, stop_gps, dataset, flag)
	global value_cache
	value_cache={}
	return(a)
	

def power_sum(c):
	ans=get_power_sum(c)
	return ans

c1=candidate();
c1.frequency=867.3;
c1.frequency=100.1;
c1.spindown=0;
c1.ra=4.904057;
c1.dec=-0.194548;
c1.iota=0;
c1.psi=0;
#c1.phi=0;

#c1.dump()


#powerflux.init("--config=search.15.config")

if 1 : powerflux.init("--dataset=random2.dst --first-bin=1953405 --side-cut=1200 -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--averaging-mode=matched --do-cutoff=0 \
		--subtract-background=0 --ks-test=1 --filter-lines=0 \
		--spindown-start=0.0 --spindown-step=2e-8 --spindown-count=1 \
		--dump-candidates=0 \
		--dump-points=0 --no-secondary-skymaps=1 \
		--sky-marks-file=all_sky_marks.txt \
		--earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
		--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--max-candidates=10000 --skymap-resolution-ratio=1 \
  		--fake-ref-time=793154935 \
  		--fake-ra=4.0904057 --fake-dec=-0.634548 --fake-freq=1085.34027777 --fake-strain=10e-24 \
  		--fake-spindown=-4.66e-10 --fake-iota=1.57 --fake-psi=0.0 \
		")
		
else :
	powerflux.init("--first-bin=1953405 --side-cut=1200 --dataset=/home/volodya/LIGO/S5/n20/dataset.H1L1.A5C.dst --config=/home/volodya/LIGO/S5/n20/matched.interactive.config")

start_center = candidate()

start_center.frequency=1085.34;
start_center.spindown=-4.5e-9;
start_center.ra=4.07;
start_center.dec=-0.62;
start_center.iota=0;
start_center.psi=0;

#start_center=c1

start_center.frequency_step=5e-5;
start_center.spindown_step=1e-12;
start_center.ra_step=0.0001;
start_center.dec_step=0.0001;
start_center.iota_step=0.01;
start_center.psi_step=0.01;




#c1.dump()
compute(start_center).dump()

datasets=powerflux.get_datasets()
#powerflux.compute_scores(c1)
#c1.dump()

#
# GUI
#

#from qtcanvas import *

class line_plot(QGraphicsView):
	def __init__(self, parent, xlab="x", ylab="y", title="", width=640, height=400):
		
		self.parent=parent
		self.app=parent.app
		
		self.curves=[]
		self.diurnal_data=[]
		self.start_gps=630720013
		self.colors=["blue", "red", "green", "black", "magenta", "cyan"]
		self.title=title
		self.xlab=xlab
		self.ylab=ylab
		
		self.x0=20
		self.y0=20
		self.width=width-2*self.x0
		self.height=height-2*self.y0
		
		self.canvas=QGraphicsScene(0, 0, self.x0*2+self.width, self.y0*2+self.height)
		QGraphicsView.__init__(self, self.canvas, parent)

		self.setFixedSize(self.x0*2+self.width+1, self.y0*2+self.height+1)
		#self.setVScrollBarMode(QScrollView.AlwaysOff)
		#self.setHScrollBarMode(QScrollView.AlwaysOff)

		self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

		self.xmax=None
		self.ymax=None
		self.xmin=None
		self.ymin=None

		self.lines=[]
		
		self.view_mode="xy"

		self.setup()
		
	def setup(self):
				
		self.canvas.setBackgroundBrush(QBrush(QColor(255, 255, 255)))
		
		
		a=self.canvas.addText(self.title)
		a.moveBy(self.x0+self.width*0.5-a.boundingRect().width()*0.5, 0)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.title_label=a
		
		a=self.canvas.addText(self.xlab)
		a.moveBy(self.x0+self.width*0.5-a.boundingRect().width()*0.5, self.y0+self.height)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.var1_label=a
		
		a=self.canvas.addText(self.ylab)
		a.moveBy(0, self.y0+self.height*0.5+a.boundingRect().width()*0.5)
		a.rotate(-90)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.var2_label=a
		
		a=self.canvas.addText(str(self.xmax))
		a.moveBy(2*self.x0+self.width-10-a.boundingRect().width(), self.height+self.y0)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.max_x_label=a

		a=self.canvas.addText(str(self.xmin))
		a.moveBy(self.x0, self.height+self.y0)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.min_x_label=a
		
		a=self.canvas.addText(str(self.ymax))
		a.rotate(-90)
		a.moveBy(0, self.y0+a.boundingRect().width())
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.max_y_label=a

		a=self.canvas.addText(str(self.ymin))
		a.rotate(-90)
		a.moveBy(0, self.height)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.min_y_label=a

	def refresh_xy(self):
		Y=[]
		X=[]
		self.points=[]
		for curve in self.curves :
			for (x,y) in curve:
				X.append(x)
				Y.append(y)
			
		self.ymax=max(Y)
		self.ymin=min(Y)
		self.xmax=max(X)
		self.xmin=min(X)
		
		delta=self.max_x_label.boundingRect().width()
		self.max_x_label.setPlainText(str(self.xmax))
		delta=delta-self.max_x_label.boundingRect().width()
		self.max_x_label.moveBy(delta, 0)
		self.min_x_label.setPlainText(str(self.xmin))
		
		delta= -self.max_y_label.boundingRect().width()
		self.max_y_label.setPlainText(str(self.ymax))
		delta=delta+self.max_y_label.boundingRect().width()
		self.max_y_label.moveBy(0, delta)
		self.min_y_label.setPlainText(str(self.ymin))
		
		for a in self.lines :
			self.canvas.removeItem(a)
		self.lines=[]
		
		dx=(self.xmax-self.xmin)
		dy=(self.ymax-self.ymin)
		
		if dx<=0 : dx=0.5
		if dy<=0 : dy=0.5
		
		color_idx=0
		for curve in self.curves:
			vec=QPainterPath()
			for (x,y) in curve:
				i=self.x0+self.width*(x-self.xmin)/dx
				j=self.y0+self.height*(self.ymax-y)/dy
				if vec.elementCount()<1:
					vec.moveTo(i,j)
				else:
					vec.lineTo(i,j)
				
			a=QGraphicsPathItem(vec)
			a.setPen(QPen(QColor(self.colors[color_idx])))
			color_idx=color_idx+1
			if color_idx>=len(self.colors) : color_idx=0
			self.canvas.addItem(a)
			self.lines.append(a)
			
	def refresh_diurnal(self):
		self.curves=[]
		for piece in self.diurnal_data:
			accum=[]
			for i in range(24):
				accum.append(0.0)
			count=copy.copy(accum)
			for (gps, x) in piece:
				hours=(gps-self.start_gps)//3600
				hours=hours % 24
				accum[hours]=accum[hours]+x
				count[hours]=count[hours]+1
			curve=[]
			for i in range(24):
				y=accum[i]/count[i]
				curve.append((i, y))
			self.curves.append(curve)
		self.refresh_xy()
			
	def refresh(self):
		if self.view_mode=="xy" : 
			self.refresh_xy()
			return
		if self.view_mode=="diurnal" : 
			self.refresh_diurnal()
			return
		print "unknown view mode "+self.view_mode

	def subtract_means(self, start=0, stop=None) :
		if stop==None : stop=len(self.curves)
		mean_curves=self.curves[start:stop]
		N=len(mean_curves)
		
		if N<1 : return
		
		for i in range(0, len(mean_curves[0])):
			sum=0
			for k in range(0, N):
				sum+=mean_curves[k][i][1]
			sum=sum/N
			for k in range(len(self.curves)):
				x=self.curves[k][i]
				self.curves[k][i]=(x[0], x[1]-sum)
		

class clickable_point(QGraphicsEllipseItem):
	def __init__(self, parent, x, y, width, height):
		QGraphicsEllipseItem.__init__(self, x, y, width, height)
		self.cand=None
		self.parent=parent
		
	def mousePressEvent(self, event):
		print "pressed"
		if(self.cand!=None):
			#self.cand.dump()
			self.cand=compute(self.cand, mode="full")
			self.parent.parent.cand_entry.setNote("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (self.cand.snr, self.cand.strain, self.cand.power_cor, self.cand.ifo_freq))
			#c=self.parent.parent.cand_entry.cand
			#setattr(c, self.parent.var1, getattr(self.cand, self.parent.var1))
			#setattr(c, self.parent.var2, getattr(self.cand, self.parent.var2))
			#self.parent.parent.cand_entry.set_controls()
			##self.update_widget.rescale()
			##self.parent.center=copy.copy(c)
			#self.parent.refresh()
			#self.parent.parent.update_plots()
	
	def mouseDoubleClickEvent(self, event):
		if(self.cand!=None):
			#self.cand.dump()
			self.parent.parent.cand_entry.setNote("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (self.cand.snr, self.cand.strain, self.cand.power_cor, self.cand.ifo_freq))
			c=self.parent.parent.cand_entry.cand
			setattr(c, self.parent.var1, getattr(self.cand, self.parent.var1))
			setattr(c, self.parent.var2, getattr(self.cand, self.parent.var2))
			self.parent.parent.cand_entry.set_controls()
			#self.update_widget.rescale()
			#self.parent.center=copy.copy(c)
			self.parent.refresh()
			self.parent.parent.update_plots()


class plot_display(QGraphicsView):
	def __init__(self, parent, var1, var2="snr"):
		
		self.parent=parent
		self.app=parent.app
		self.var1=var1
		self.var2=var2
		
		self.dot_half_count=20
		
		
		self.x0=20
		self.y0=20
		self.width=100*4
		self.height=100*4
		
		self.canvas=QGraphicsScene(0, 0, self.x0*2+self.width, self.y0*2+self.height)
		QGraphicsView.__init__(self, self.canvas, parent)

		self.setFixedSize(self.x0*2+self.width+1, self.y0*2+self.height+1)
		#self.setVScrollBarMode(QScrollView.AlwaysOff)
		#self.setHScrollBarMode(QScrollView.AlwaysOff)

		self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

		self.xmax=None
		self.ymax=None
		self.xmin=None
		self.ymin=None

		self.dots=[]

		self.setup()
		
	def setup(self):
		
		#self.layout = QGridLayout(self, 1, 1)
		
		#self.canvas_view=QCanvasView(self.canvas, self)
		#self.canvas_view.setFixedSize(size, size)
		#self.canvas_view.setVScrollBarMode(QScrollView.AlwaysOff)
		#self.canvas_view.setHScrollBarMode(QScrollView.AlwaysOff)
		
		#self.canvas_view.setCanvas(QCanvas())
		
		#self.layout.addWidget(self.canvas_view, 0, 0)
		
		self.canvas.setBackgroundBrush(QBrush(QColor(255, 255, 255)))
		
		#self.canvas.setAdvancePeriod(30)
		
		a=self.canvas.addText(self.var1)
		a.moveBy(self.x0+self.width*0.5, self.y0+self.height)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.var1_label=a
		
		a=self.canvas.addText(self.var2)
		a.moveBy(0, self.y0+self.height*0.5)
		a.rotate(-90)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.var2_label=a
		
		a=self.canvas.addText(str(self.xmax))
		a.moveBy(2*self.x0+self.width-10-a.boundingRect().width(), self.height+self.y0)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.max_x_label=a

		a=self.canvas.addText(str(self.xmin))
		a.moveBy(self.x0, self.height+self.y0)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.min_x_label=a
		
		a=self.canvas.addText(str(self.ymax))
		a.rotate(-90)
		a.moveBy(0, self.y0+a.boundingRect().width())
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.max_y_label=a

		a=self.canvas.addText(str(self.ymin))
		a.rotate(-90)
		a.moveBy(0, self.height-self.y0)
		a.setDefaultTextColor(QColor(0,0,0))
		a.show()
		self.min_y_label=a

	def refresh(self):
		n=self.dot_half_count
		c=copy.copy(self.center)
		Y=[]
		X=[]
		self.points=[]
		for i in range(-n, n+1) :
			setattr(c, self.var1, getattr(self.center, self.var1)+i*getattr(self.center, self.var1+"_step"))
			X.append(getattr(c, self.var1))
			a=compute(c)
			y=getattr(a, self.var2)
			Y.append(y)
			self.points.append(a)
			
		self.ymax=max(Y)
		self.ymin=min(Y)
		self.xmax=max(X)
		self.xmin=min(X)
		
		delta=self.max_x_label.boundingRect().width()
		self.max_x_label.setPlainText(str(self.xmax))
		delta=delta-self.max_x_label.boundingRect().width()
		self.max_x_label.moveBy(delta, 0)
		self.min_x_label.setPlainText(str(self.xmin))
		
		delta= -self.max_y_label.boundingRect().width()
		self.max_y_label.setPlainText(str(self.ymax))
		delta=delta+self.max_y_label.boundingRect().width()
		self.max_y_label.moveBy(0, delta)
		self.min_y_label.setPlainText(str(self.ymin))
		
		for a in self.dots :
			self.canvas.removeItem(a)
		self.dots=[]
		
		dx=(self.xmax-self.xmin)
		dy=(self.ymax-self.ymin)
		
		if dx<=0 : dx=0.5
		if dy<=0 : dy=0.5
		
		for (x,y, c) in map(None, X, Y, self.points) :
			i=self.x0+self.width*(x-self.xmin)/dx
			j=self.y0+self.height*(self.ymax-y)/dy
			a=clickable_point(self, i, j, 8, 8)
			if y>=self.ymax :
				a.setBrush(QBrush(QColor(255,0,0)))
			else :
				a.setBrush(QBrush(QColor(100,100,255)))
			a.cand=c
			self.canvas.addItem(a)
			self.dots.append(a)
		app.processEvents()
		
	def rescale(self):
		pass
		
	def best_candidate(self):
		for c in self.points :
			if(getattr(c, self.var2)>=self.ymax) : return c
		return None
		
	def update(self, center):
		self.center=center
		self.refresh()


class clickable_dot(QGraphicsRectItem):
	def __init__(self, parent):
		QGraphicsRectItem.__init__(self)
		self.cand=None
		self.parent=parent
		self.update_widget=parent.update_widget
		
	def mousePressEvent(self, event):
		print "pressed"
		if(self.cand!=None):
			#self.cand.dump()
			self.cand=compute(self.cand, mode="full")
			self.update_widget.parent.cand_entry.setNote("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (self.cand.snr, self.cand.strain, self.cand.power_cor, self.cand.ifo_freq))
			c=self.update_widget.parent.cand_entry.cand
			setattr(c, self.update_widget.var1, getattr(self.cand, self.update_widget.var1))
			setattr(c, self.update_widget.var2, getattr(self.cand, self.update_widget.var2))
			self.update_widget.parent.cand_entry.set_controls()
			self.parent.subdivide()
			self.update_widget.rescale()
			self.update_widget.refresh()
	
	def mouseDoubleClickEvent(self, event):
		if(self.cand!=None):
			c=copy.copy(self.cand)
			self.update_widget.parent.cand_entry.setNote("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (self.cand.snr, self.cand.strain, self.cand.power_cor, self.cand.ifo_freq))
			c1=self.update_widget.parent.cand_entry.cand
			setattr(c1, self.update_widget.var1, getattr(self.cand, self.update_widget.var1))
			setattr(c1, self.update_widget.var2, getattr(self.cand, self.update_widget.var2))
			self.update_widget.parent.cand_entry.set_controls()
			self.update_widget.update(c)
			self.update_widget.parent.update(c)

class divisible_dot:
	def __init__(self, update_widget, center, x, y, size, var1, var2, var3="snr"):
		self.center=compute(center)
		self.update_widget=update_widget
		self.size=size
		self.x=x
		self.y=y
		self.var1=var1
		self.var2=var2
		self.var3=var3
		
		self.step1=getattr(center, self.var1+"_step")
		self.step2=getattr(center, self.var2+"_step")

		self.dots=None
				
		a=clickable_dot(self)
		self.mydot=a
		a.setRect(self.x, self.y, self.size, self.size)
		a.cand=self.center
		self.update_widget.canvas.addItem(a)
		self.update()
		
	def __del__(self):
		self.mydot.hide()
		self.update_widget.canvas.removeItem(self.mydot)
		
	def hide(self):
		self.mydot.hide()
		if(self.dots!=None) :
			for b in self.dots: b.hide()
	
		
	def update(self):
		
		if(self.dots==None):
			a=self.mydot
			z=255*(getattr(self.center, self.var3)-self.update_widget.min_z)/(self.update_widget.max_z-self.update_widget.min_z)
			if(z<0):z=0
			if(z>255):z=255
			z=255-z
			if(self.center.strain> -0.99):
				a.setBrush(QBrush(QColor.fromHsv(z, 255, 255)))
				a.setPen(QPen(QColor.fromHsv(z, 255, 255)))
			else:
				a.setBrush(QBrush(QColor(0, 0, 0)))
				a.setPen(QPen(QColor(0, 0, 0)))
			a.show()
		else:
			for i in range(0, len(self.dots)):
				self.dots[i].update()
		
	def max(self):
		if(self.dots==None):
			return getattr(self.center, self.var3)
		else:
			return max([a.max() for a in self.dots])
		
	def min(self):
		if(self.dots==None):
			return getattr(self.center, self.var3)
		else:
			return min([a.min() for a in self.dots])

	def auto_subdivide1(self, level, half_n=2):
		if level < 0 : level=self.max()
		a=level
		for b in self.dots:
			if b.max()>=level : 
				b.subdivide(half_n=half_n, auto_level=None)
				x=b.min()
				if x<a : 
					a=x
				
		
		
		for b in self.dots:
			if b.max() > a : b.subdivide(half_n=half_n, auto_level=None)

	def auto_subdivide(self, level, half_n=2):
		a=[b.max() for b in self.dots]
		a.sort()
		a.reverse()
		a=a[3]
		
		for b in self.dots:
			if b.max() > a : b.subdivide(half_n=half_n, auto_level=None)

	def subdivide(self, half_n=2, auto_level=None):
		if(self.dots!=None) :
			if(auto_level!=None) : self.auto_subdivide(half_n=half_n, level=auto_level)
			return
			
		
		self.half_n=half_n
		n=half_n
		s=self.size/(2*n+1)
		if s<2 : return
		
		c=copy.copy(self.center)
		setattr(c, self.var1+"_step", self.step1/(2*n+1))
		setattr(c, self.var2+"_step", self.step2/(2*n+1))
		self.dots=[]
		for i in range(-n, n+1) :
			setattr(c, self.var1, getattr(self.center, self.var1)+i*self.step1)
			for j in range(-n, n+1) :
				setattr(c, self.var2, getattr(self.center, self.var2)+j*self.step2)
				a=divisible_dot(self.update_widget, c, self.x+(i+n)*s, self.y+(j+n)*s, s, self.var1, self.var2)
				self.dots.append(a)
		self.mydot.hide()
		
		self.update_widget.app.processEvents()

		if(auto_level!=None) : self.auto_subdivide(half_n=half_n, level=auto_level)
		
	def best_candidate(self):
		if(self.dots==None):
			return self.center
		c=self.center
		for b in self.dots :
			x=b.max()
			if x> getattr(c, self.var3) : c=b.best_candidate()
		return c


class map_display(QGraphicsView):
	def __init__(self, parent, var1, var2, var3="snr"):
		
		self.parent=parent
		self.app=parent.app
		self.var1=var1
		self.var2=var2
		self.var3=var3
		
		self.dot_size=20
		self.dot_half_count=10
		
		self.max_z=10
		self.min_z=0
		
		#self.size=(2*self.dot_half_count+1)*self.dot_size+2*self.dot_size
		self.size=2*self.dot_size + 4*100
		self.canvas=QGraphicsScene(0, 0, self.size, self.size)
		QGraphicsView.__init__(self, self.canvas, parent)

		self.setFixedSize(self.size+1, self.size+1)
		#self.setVScrollBarMode(QScrollView.AlwaysOff)
		#self.setHScrollBarMode(QScrollView.AlwaysOff)

		self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

		self.setup()
		
	def setup(self):
		s=self.dot_size
		n=self.dot_half_count
		
		#self.layout = QGridLayout(self, 1, 1)
		
		#self.canvas_view=QCanvasView(self.canvas, self)
		#self.canvas_view.setFixedSize(size, size)
		#self.canvas_view.setVScrollBarMode(QScrollView.AlwaysOff)
		#self.canvas_view.setHScrollBarMode(QScrollView.AlwaysOff)
		
		#self.canvas_view.setCanvas(QCanvas())
		
		#self.layout.addWidget(self.canvas_view, 0, 0)
		
		self.canvas.setBackgroundBrush(QBrush(QColor(0, 0, 0)))
		
		#self.canvas.setAdvancePeriod(30)
		
		a=self.canvas.addText(self.var1)
		a.moveBy(self.size*0.5, self.dot_size*0)
		a.setDefaultTextColor(QColor(255,255,255))
		a.show()
		self.var1_label=a
		
		a=self.canvas.addText(self.var2)
		a.moveBy(self.dot_size*0, self.size*0.5)
		a.rotate(-90)
		a.setDefaultTextColor(QColor(255,255,255))
		a.show()
		self.var2_label=a
		
		a=self.canvas.addText(str(self.max_z))
		a.moveBy(self.size-10-a.boundingRect().width(), self.size-self.dot_size*1.0)
		a.setDefaultTextColor(QColor(255,255,255))
		a.show()
		self.max_z_label=a

		a=self.canvas.addText(str(self.min_z))
		a.moveBy(10, self.size-self.dot_size*1.0)
		a.setDefaultTextColor(QColor(255,255,255))
		a.show()
		self.min_z_label=a

		self.dot=None

		
	def refresh(self):
		self.dot.update()
		
	def rescale(self):
		self.max_z=self.dot.max()
		self.min_z=self.dot.min()
		
		if self.max_z<=self.min_z :
			a=0.01*self.max_z
			if(a<=0) : a=0.01
			self.min_z=self.max_z-a
		
		delta=self.max_z_label.boundingRect().width()
		self.max_z_label.setPlainText(str(self.max_z))
		delta=delta-self.max_z_label.boundingRect().width()
		self.max_z_label.moveBy(delta, 0)
		self.min_z_label.setPlainText(str(self.min_z))
		
	def best_candidate(self):
		if self.dot!=None : return self.dot.best_candidate()
		return None
		
	def update(self, center):
		
		s=self.dot_size
		n=self.dot_half_count
		self.center=center
		
		if(self.dot!=None) :
			self.dot.hide()
			del self.dot
		self.dot=divisible_dot(self, self.center, s, s, 4*100, self.var1, self.var2, self.var3)
		
		self.dot.subdivide(auto_level=-1)
		self.rescale()
		self.dot.update()
		
candidate_items=["frequency", "spindown", "ra", "dec", "iota", "psi"]
		
class candidate_display(QWidget):
	def __init__(self, parent):
		self.parent=parent
		self.cand=None
		QWidget.__init__(self, parent)
		
		layout=QGridLayout(self)
		self.layout=layout
	
		self.display=QTableWidget(2, len(candidate_items))
		
		for i in range(len(candidate_items)):
			self.display.setHorizontalHeaderItem(i, QTableWidgetItem(candidate_items[i]))
		self.display.setVerticalHeaderItem(0, QTableWidgetItem("value"))
		self.display.setVerticalHeaderItem(1, QTableWidgetItem("step"))
		
		layout.addWidget(self.display, 0, 0)
	
		self.note=QLabel("", self)
		layout.addWidget(self.note, 0, 1)

		self.controls=QWidget(self)
		layout.addWidget(self.controls, 0, 2)
	
		layout2=QGridLayout(self.controls)
	
		self.compute = QPushButton("Compute", self)
		layout2.addWidget(self.compute, 0, 0)
		app.connect(self.compute, SIGNAL("clicked()"), self.recompute)
		
		self.recompute_all = QPushButton("Recompute all", self)
		layout2.addWidget(self.recompute_all, 1, 0)
		app.connect(self.recompute_all, SIGNAL("clicked()"), self.parent.recompute)

	def recompute(self):
		self.get_controls()
		self.cand=compute(self.cand, mode="full")
		c=self.cand
		self.setNote("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (c.snr, c.strain, c.power_cor, c.ifo_freq))

	def set_controls(self):
		if(self.cand==None): return
		
		for i in range(len(candidate_items)):
			value=getattr(self.cand, candidate_items[i], "")
			self.display.setItem(0, i, QTableWidgetItem(str(value)))

			value=getattr(self.cand, candidate_items[i]+"_step", "")
			self.display.setItem(1, i, QTableWidgetItem(str(value)))

	def get_controls(self):
		
		for i in range(len(candidate_items)):
			value=float(self.display.item(0, i).text())
			setattr(self.cand, candidate_items[i], value)

			value=float(self.display.item(1, i).text())
			setattr(self.cand, candidate_items[i]+"_step", value)
			
	def set_coords(self, frequency=None, spindown=None, ra=None, dec=None, iota=None, psi=None):
		if frequency!=None : cand.frequency=frequency
		if spindown!=None : cand.spindown=spindown
		if ra!=None : cand.ra=ra
		if dec!=None : cand.dec=dec
		if iota!=None : cand.iota=iota
		self.set_controls()
		
	def setNote(self, string):
		self.note.setText(string)

class powerflux_main_window(QMainWindow):
    """An application called astrocalc."""

    def __init__(self, app):
        QMainWindow.__init__(self, None)
        #self.initIcons()
	self.app=app
        self.setup()
        #self.initPrinter()
        #self.initToolBar()
        #self.initMenu()
        self.initMainWidget()
        #self.setCaption(self.appTitle)

    def setup(self):
        self.appTitle = "PowerFlux - study CW signals"
	
	self.center = start_center
	
    def set_controls(self):
	self.cand_entry.cand=self.center
	self.cand_entry.set_controls()

    def get_controls(self):
	self.cand_entry.get_controls()
	self.center=self.cand_entry.cand

    def initMainWidget(self):
	    
	self.main=QWidget(self)
	self.setCentralWidget(self.main)
	
	layout0 = QGridLayout(self.main)
		
	self.cand_entry=candidate_display(self)
	#self.controls.layout2.addWidget(self.cand_entry, 0, 0, 1, 3)
	layout0.addWidget(self.cand_entry, 1, 0)
	
	#self.f=QWidget(self)
	self.tabs=QTabWidget(self)
	#self.setCentralWidget(self.tabs)
	layout0.addWidget(self.tabs, 0, 0)
	
	app.connect(self.tabs, SIGNAL("currentChanged(int)"), self.tab_switch)
	
	f=QWidget()
	self.tabs.addTab(f, "Maps")
	#QGridLayout(self.tabs.widget(0)).addWidget(f)
	self.layout1=QGridLayout(f)
	
	self.controls = QWidget(self)
	self.layout1.addWidget(self.controls, 1, 2)
	
	
	self.controls.layout2 = QGridLayout(self.controls)
	
	#self.controls.recompute = QPushButton("Recompute all", self.controls)
	#self.controls.layout2.addWidget(self.controls.recompute, 8, 0, 1, 3)
	#app.connect(self.controls.recompute, SIGNAL("clicked()"), self.recompute)

	self.controls.best_candidate = QPushButton("Best candidate", self.controls)
	self.controls.layout2.addWidget(self.controls.best_candidate, 9, 0, 1, 3)
	app.connect(self.controls.best_candidate, SIGNAL("clicked()"), self.best_candidate)
	
	self.set_controls()

	self.maps=[]
	
	f=map_display(self, "frequency", "spindown")
	self.layout1.addWidget(f, 0, 0)
	self.maps.append(f)
	
	f=map_display(self, "ra", "dec")
	self.layout1.addWidget(f, 0, 1)
	self.maps.append(f)

	f=map_display(self, "spindown", "dec")
	self.layout1.addWidget(f, 1, 0)
	self.maps.append(f)

	f=map_display(self, "spindown", "ra")
	self.layout1.addWidget(f, 1, 1)
	self.maps.append(f)

	f=map_display(self, "iota", "psi")
	self.layout1.addWidget(f, 0, 2)
	self.maps.append(f)


	f=QWidget()
	self.tabs.addTab(f, "Graphs")
	layout=QGridLayout(f)

	self.plots=[]
	
	f=plot_display(self, "frequency")
	layout.addWidget(f, 0, 0)
	self.plots.append(f)

	f=plot_display(self, "ra")
	layout.addWidget(f, 0, 1)
	self.plots.append(f)

	f=plot_display(self, "dec")
	layout.addWidget(f, 0, 2)
	self.plots.append(f)

	f=plot_display(self, "spindown")
	layout.addWidget(f, 1, 0)
	self.plots.append(f)

	f=plot_display(self, "iota")
	layout.addWidget(f, 1, 1)
	self.plots.append(f)

	f=plot_display(self, "psi")
	layout.addWidget(f, 1, 2)
	self.plots.append(f)

	f=QWidget()
	self.tabs.addTab(f, "Datasets")
	layout=QGridLayout(f)

	self.dataset_plots=[]
	
	f=line_plot(self, "frequency", "average power", "new_weighted_mean", width=640, height=440)
	layout.addWidget(f, 0, 0)
	self.dataset_plots.append(f)
	
	curves=[]
	for dataset in datasets:
		curve=[]
		step=1.0/dataset['coherence_time']
		frequency=dataset['first_bin']/dataset['coherence_time']
		for power in dataset['new_weighted_means']:
			curve.append((frequency, power))
			frequency=frequency+step
		curves.append(curve)		
	f.curves=curves
	
	f=line_plot(self, "frequency", "log(power)", "FMedians", width=640, height=440)
	layout.addWidget(f, 1, 0)
	self.dataset_plots.append(f)
	
	curves=[]
	for dataset in datasets:
		curve=[]
		step=1.0/dataset['coherence_time']
		frequency=dataset['first_bin']/dataset['coherence_time']
		for power in dataset['FMedians']:
			curve.append((frequency, power))
			frequency=frequency+step
		curves.append(curve)		
	f.curves=curves

	f=line_plot(self, "hours", "log(power)", "TMedians", width=640, height=440)
	layout.addWidget(f, 0, 1)
	self.dataset_plots.append(f)
	
	curves=[]
	for dataset in datasets:
		curve=[]
		
		start=dataset['gps'][0]
		for (gps,power) in map(None, dataset['gps'], dataset['TMedians']):
			time=(gps-start)/3600.0
			curve.append((time, power))
		curves.append(curve)		
	f.curves=curves
	
	rows=["name", "detector", "sft_count", "weight", "nbins", "coherence_time", "gps_start", "gps_stop"]
	
	f=QTableWidget(len(rows), len(datasets)+1)
	layout.addWidget(f, 1, 1)
	
	
	for i in range(0, len(rows)):
		f.setItem(i, 0, QTableWidgetItem(rows[i]))
		
	for j in range(0, len(datasets)):
		dataset=datasets[j]
		for i in range(0, len(rows)):
			f.setItem(i, j+1, QTableWidgetItem(str(dataset[rows[i]])))
		
	f.resizeColumnsToContents()
		
	f=QWidget()
	self.tabs.addTab(f, "SFT info")
	layout=QGridLayout(f)
	
	rows=0
	for dataset in datasets:
		rows+=len(dataset['gps'])
	
	self.row_info_columns=["gps", "dataset", "detector", "TMedian", "sft_veto", "det_vel[0]", "det_vel[1]", "det_vel[2]", "weight", "power"]
	self.row_info=QTableWidget(rows, len(self.row_info_columns))
	layout.addWidget(self.row_info, 0, 0)
	
	
	f=QWidget()
	self.tabs.addTab(f, "Power accumulation")
	layout=QGridLayout(f)

	self.power_accum_plots=[]
	
	self.power_evolution=line_plot(self, "Nsft", "Power", "Power evolution by SFT number", width=640, height=440)
	layout.addWidget(self.power_evolution, 0, 0)
	self.power_accum_plots.append(self.power_evolution)

	self.power_evolution2=line_plot(self, "weight", "Power", "Power evolution by weight", width=640, height=440)
	layout.addWidget(self.power_evolution2, 1, 0)
	self.power_accum_plots.append(self.power_evolution2)

	self.avg_power_evolution=line_plot(self, "Nsft", "Power", "Average power evolution by SFT number", width=640, height=440)
	layout.addWidget(self.avg_power_evolution, 0, 1)
	self.power_accum_plots.append(self.avg_power_evolution)

	f=QWidget()
	self.tabs.addTab(f, "Diurnal variation")
	layout=QGridLayout(f)

	self.diurnal_plots=[]
	
	self.diurnal_tmedians=line_plot(self, "Hour", "Power", "Diurnal variation of SFT noise level", width=640, height=440)
	layout.addWidget(self.diurnal_tmedians, 0, 0)
	self.diurnal_plots.append(self.diurnal_tmedians)

	self.diurnal_weight=line_plot(self, "Hour", "Weight", "Diurnal variation of SFT weight", width=640, height=440)
	layout.addWidget(self.diurnal_weight, 1, 0)
	self.diurnal_plots.append(self.diurnal_weight)

	self.diurnal_power=line_plot(self, "Hour", "Power", "Diurnal variation of power", width=640, height=440)
	layout.addWidget(self.diurnal_power, 0, 1)
	self.diurnal_plots.append(self.diurnal_power)
	
	self.diurnal_summand=line_plot(self, "Hour", "Power", "Diurnal variation of power sum components", width=640, height=440)
	layout.addWidget(self.diurnal_summand, 1, 1)
	self.diurnal_plots.append(self.diurnal_summand)

	for plot in self.diurnal_plots:
		plot.view_mode="diurnal"

	f=QWidget()
	self.tabs.addTab(f, "Freehand")
	layout=QGridLayout(f)
	
	self.freehand_output=QTextEdit()
	layout.addWidget(self.freehand_output, 0, 0)
	#self.freehand_output.readOnly=1
	
	self.freehand_input=QLineEdit()
	layout.addWidget(self.freehand_input, 1, 0)
	app.connect(self.freehand_input, SIGNAL("returnPressed()"), self.execute_freehand_input)

	#self.update(self.center)
	
    def write(self, string, color="black"):
	#cursor=self.freehand_output.currentCursor()
	self.freehand_output.moveCursor(QTextCursor.End, QTextCursor.MoveAnchor)
	self.freehand_output.setTextColor(QColor(color))
	self.freehand_output.insertPlainText(string)
	
    def execute_freehand_input(self):
	cmd=self.freehand_input.text()
	self.write(cmd+"\n", color="blue")
	print cmd
	stdout=sys.stdout
	sys.stdout=self
	code=compile(str(cmd), "<freehand>", 'single')
	eval(code)
	sys.stdout=stdout
	self.write("\n")
	
    def tab_switch(self, index):
	print "switched"
	self.recompute()
	#print index
	#self.tabs.widget(index).layout().addWidget(self.controls, 1, 2)

    def best_candidate(self):
	    c=None
	    for w in self.maps :
		    c2=w.best_candidate()
		    if c2 == None : continue
		    if c==None :
			    c=c2
			    continue
		    if c2.snr>c.snr : c=c2
	    if c==None : return
	    self.cand_entry.cand.copy_coords(c)
	    self.cand_entry.set_controls()
	    self.cand_entry.recompute()

    def recompute(self):
	    self.cand_entry.get_controls()
	    
	    index=self.tabs.tabText(self.tabs.currentIndex())
	    print index
	    if index=="Maps" :
			self.update_maps()
			return
	    if index=="Graphs":
			self.update_plots()
			return
	    if index=="Datasets":
			self.update_datasets()
			return
	    if index=="SFT info":
			self.update_row_info()
			return
	    if index=="Power accumulation":
			self.update_power_accumulation()
			return
	    if index=="Diurnal variation":
			self.update_diurnal_variation()
			return

    def update_maps(self, center=None):
	    if center == None : center=self.cand_entry.cand
	    for plot in self.maps :
		    plot.update(center)
	    
    def update_plots(self, center=None):
	    if center == None : center=self.cand_entry.cand
	    for plot in self.plots :
		    plot.update(center)
	
    def update_datasets(self, center=None):
	    if center == None : center=self.cand_entry.cand
  	    for plot in self.dataset_plots:
		  plot.refresh()

    def update_power_accumulation(self, center=None):
	    if center == None : center=self.cand_entry.cand
	    
	    psums=[power_sum(center)]
	    c=copy.copy(center)
	    for i in range(3, 10) :
		    c.frequency=center.frequency+i*0.01
		    psums.append(power_sum(c))
		    c.frequency=center.frequency-i*0.01
		    psums.append(power_sum(c))
	    
	    self.power_evolution.curves=[]
	    self.power_evolution2.curves=[]
	    self.avg_power_evolution.curves=[]
	    
	    for psum in psums :
		power_curve=[]
		power_curve2=[]
		avg_power_curve=[]
		weight=0
		sum=0
		x=0
		for (w,p) in psum :
			weight=weight+w
			sum=sum+p*w
			x=x+1
			power=sum/weight
			power_curve.append((x, sum))
			power_curve2.append((weight, sum))
			avg_power_curve.append((x, power))
			
		self.power_evolution.curves.append(power_curve)

		self.power_evolution2.curves.append(power_curve2)
		
		self.avg_power_evolution.curves.append(avg_power_curve)
	    
  	    self.power_evolution.subtract_means(start=1)
     	    self.power_evolution2.subtract_means(start=1)

  	    for plot in self.power_accum_plots:
		  plot.refresh()
		  
    def update_row_info(self, center=None):
	if center == None : center=self.cand_entry.cand
	
	for i in range(len(self.row_info_columns)):
		self.row_info.setHorizontalHeaderItem(i, QTableWidgetItem(self.row_info_columns[i]))
		
	row=0
	for dataset in datasets:
		gps=dataset['gps']
		dname=dataset['name']
		detector=dataset['detector']
		TMedian=dataset['TMedians']
		det_vel=dataset['detector_velocity']
		veto=dataset['sft_veto']
		
		for i in range(len(gps)):
			self.row_info.setItem(row, 0, QTableWidgetItem(str(gps[i])))
			self.row_info.setItem(row, 1, QTableWidgetItem(dname))
			self.row_info.setItem(row, 2, QTableWidgetItem(detector))
			self.row_info.setItem(row, 3, QTableWidgetItem(str(TMedian[i])))
			self.row_info.setItem(row, 4, QTableWidgetItem(str(veto[i])))
			
			self.row_info.setItem(row, 5, QTableWidgetItem(str(det_vel[i][0])))
			self.row_info.setItem(row, 6, QTableWidgetItem(str(det_vel[i][1])))
			self.row_info.setItem(row, 7, QTableWidgetItem(str(det_vel[i][2])))
			row+=1
		
	psum=power_sum(center)
	
	for i in range(len(psum)):
		x=psum[i]
		self.row_info.setItem(i, 8, QTableWidgetItem(str(x[0])))
		self.row_info.setItem(i, 9, QTableWidgetItem(str(x[1])))
	
	self.row_info.resizeColumnsToContents()
    
    def update_diurnal_variation(self, center=None):
	if center == None : center=self.cand_entry.cand

	
	self.diurnal_tmedians.diurnal_data=[]
	for dataset in datasets :
		curve=map(None, dataset['gps'], dataset['TMedians'])
		self.diurnal_tmedians.diurnal_data.append(curve)
		
    	psum=power_sum(center)
	weight=map(lambda x: x[0], psum)
	power=map(lambda x: x[1], psum)
	summand=map(lambda x: x[0]*x[1], psum)
	
	self.diurnal_weight.diurnal_data=[]
	self.diurnal_summand.diurnal_data=[]
	self.diurnal_power.diurnal_data=[]
	
	start=0
	for dataset in datasets :
		gps=dataset['gps']
		count=len(gps)
		
		curve=map(None, gps, weight[start:start+count])
		self.diurnal_weight.diurnal_data.append(curve)
	
		curve=map(None, gps, summand[start:start+count])
		self.diurnal_summand.diurnal_data.append(curve)
		
		curve=map(None, gps, power[start:start+count])
		self.diurnal_power.diurnal_data.append(curve)
		
		start+=count
	
	
	for plot in self.diurnal_plots:
		plot.refresh()


app=None
	

def main(args):
    global app
    app=QApplication(args)
    mainWindow = powerflux_main_window(app)
    mainWindow.show()
    app.connect(app, SIGNAL("lastWindowClosed()"), app, SLOT("quit()"))
    #mainWindow.update(mainWindow.center)
    app.exec_()



if __name__ == "__main__":
    main(sys.argv)

