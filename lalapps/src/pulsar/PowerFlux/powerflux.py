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

def compute(c):
	ans=value_cache.get(c.coords())
	if(ans==None):
		ans=copy.copy(c)
		compute_scores(ans)
		value_cache[c.coords()]=ans
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

#powerflux.init("--dataset=random.dst -f 180002 -n 501 \
		#--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		#--averaging-mode=matched --do-cutoff=1 --lock-file=test.lock \
		#--subtract-background=0 --ks-test=1 --filter-lines=0 \
		#--spindown-start=0.0 --spindown-step=2e-8 --spindown-count=1 \
		#--dump-candidates=0 \
		#--dump-points=0 --no-secondary-skymaps=1 \
		#--sky-marks-file=all_sky_marks.txt \
		#--earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
		#--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		#--max-candidates=10000 --skymap-resolution-ratio=1")
		
powerflux.init("--first-bin=1953405 --side-cut=2400 --dataset=/home/volodya/LIGO/S5/n20/dataset.H1L1.A5.dst --config=/home/volodya/LIGO/S5/n20/matched.interactive.config")

start_center = candidate()

start_center.frequency=1085.34;
start_center.spindown=-4.5e-9;
start_center.ra=4.07;
start_center.dec=-0.62;
start_center.iota=0;
start_center.psi=0;

#start_center=c1

start_center.frequency_step=1e-4;
start_center.spindown_step=5e-12;
start_center.ra_step=0.001;
start_center.dec_step=0.001;
start_center.iota_step=0.1;
start_center.psi_step=0.1;


#c1.dump()
compute(c1).dump()
#powerflux.compute_scores(c1)
#c1.dump()

#
# GUI
#

#from qtcanvas import *

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
			self.update_widget.parent.controls.display.setText("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (self.cand.snr, self.cand.strain, self.cand.power_cor, self.cand.ifo_freq))
			setattr(self.update_widget.parent.center, self.update_widget.var1, getattr(self.cand, self.update_widget.var1))
			setattr(self.update_widget.parent.center, self.update_widget.var2, getattr(self.cand, self.update_widget.var2))
			self.update_widget.parent.set_controls()
			self.parent.subdivide()
			self.update_widget.rescale()
	
	def mouseDoubleClickEvent(self, event):
		if(self.cand!=None):
			c=copy.copy(self.cand)
			self.update_widget.parent.controls.display.setText("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (self.cand.snr, self.cand.strain, self.cand.power_cor, self.cand.ifo_freq))
			setattr(self.update_widget.parent.center, self.update_widget.var1, getattr(self.cand, self.update_widget.var1))
			setattr(self.update_widget.parent.center, self.update_widget.var2, getattr(self.cand, self.update_widget.var2))
			self.update_widget.parent.set_controls()
			self.update_widget.update(c)
			self.update_widget.parent.update(c)

class divisible_dot:
	def __init__(self, update_widget, center, x, y, size, var1, var2):
		self.center=compute(center)
		self.update_widget=update_widget
		self.size=size
		self.x=x
		self.y=y
		self.var1=var1
		self.var2=var2
		
		self.step1=getattr(center, self.var1+"_step")
		self.step2=getattr(center, self.var2+"_step")

		self.dots=None
				
		a=clickable_dot(self)
		self.mydot=a
		a.setRect(self.x, self.y, self.size, self.size)
		a.cand=self.center
		self.update_widget.canvas.addItem(a)
		self.update()
		
	def update(self):
		
		if(self.dots==None):
			a=self.mydot
			z=255*(getattr(self.center, "snr")-self.update_widget.min_z)/(self.update_widget.max_z-self.update_widget.min_z)
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
			return self.center.snr
		else:
			return max([a.max() for a in self.dots])
		
	def min(self):
		if(self.dots==None):
			return self.center.snr
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
			if x>c.snr : c=b.best_candidate()
		return c

class two_display(QGraphicsView):
	def __init__(self, parent, var1, var2):
		
		self.parent=parent
		self.app=parent.app
		self.var1=var1
		self.var2=var2
		
		self.dot_size=20
		self.dot_half_count=10
		
		self.max_z=10
		self.min_z=0
		
		#self.size=(2*self.dot_half_count+1)*self.dot_size+2*self.dot_size
		self.size=2*self.dot_size + 4*125
		self.canvas=QGraphicsScene(0, 0, self.size, self.size)
		QGraphicsView.__init__(self, self.canvas, parent)

		self.setFixedSize(self.size, self.size)
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

		return
		self.dots=[]
		
		
		for i in range(-n, n) :
			for j in range(-n, n) :
				a=clickable_dot(self)
				a.setRect((i+11)*s, (j+11)*s, s, s)
				self.canvas.addItem(a)
				a.setBrush(QBrush(QColor(0, 0, 0)))
				a.show()
				self.dots.append(a)
		
		#a=QCanvasRectangle(0, 0, 10, 10, self.canvas)
		#a.setBrush(QBrush(QColor(255, 255, 255)))
		#a.setPen(QPen(QColor(0.5, 0.5, 0.5)))
		#a.show()
		
	def refresh(self):
		pass
		
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
		del self.dot
		
		s=self.dot_size
		n=self.dot_half_count
		self.center=center
		
		self.dot=divisible_dot(self, self.center, s, s, 4*125, self.var1, self.var2)
		
		self.dot.subdivide(auto_level=-1)
		self.rescale()
		self.dot.update()
		
		return
	
		c=copy.copy(center)
		s=self.dot_size
		n=self.dot_half_count
		step1=getattr(center, self.var1+"_step")
		step2=getattr(center, self.var2+"_step")
		
		print self.var1+" "+self.var2
		self.points=[]
		k=0
		for i in range(-n, n) :
			setattr(c, self.var1, getattr(center, self.var1)+i*step1)
			for j in range(-n, n) :
				setattr(c, self.var2, getattr(center, self.var2)+j*step2)
				# magic happens
				#c.snr=(2*n+i+j)/4.0
				#c.strain=1
				a=self.dots[k]
				a.cand=compute(c)
				self.points.append(a.cand)
				#print c.frequency+" "+c.spindown
				z=255*(getattr(a.cand, "snr")-self.min_z)/(self.max_z-self.min_z)
				if(z<0):z=0
				if(z>255):z=255
				z=255-z
				if(a.cand.strain> -0.99):
					a.setBrush(QBrush(QColor.fromHsv(z, 255, 255)))
				else:
					a.setBrush(QBrush(QColor(0, 0, 0)))
				k=k+1
			self.app.processEvents()
		
		self.min_z=self.points[0].snr
		self.max_z=self.points[0].snr
		for k in range(2, len(self.points)):
			x=self.points[k].snr
			if(x>self.max_z) : self.max_z=x
			if(x<self.min_z) : self.min_z=x
		if(self.max_z<=self.min_z) : self.max_z=self.min_z+1
		delta=self.max_z_label.boundingRect().width()
		self.max_z_label.setPlainText(str(self.max_z))
		delta=delta-self.max_z_label.boundingRect().width()
		self.max_z_label.moveBy(delta, 0)
		self.min_z_label.setPlainText(str(self.min_z))
		k=0
		for i in range(-n, n) :
			for j in range(-n, n) :
				c=self.points[k]
				a=self.dots[k]
				z=255*(getattr(c, "snr")-self.min_z)/(self.max_z-self.min_z)
				if(z<0):z=0
				if(z>255):z=255
				z=255-z
				if(c.strain> -0.99):
					a.setBrush(QBrush(QColor.fromHsv(z, 255, 255)))
				else:
					a.setBrush(QBrush(QColor(0, 0, 0)))
				k=k+1
		self.app.processEvents()

class powerflux(QMainWindow):
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
	self.controls.entry_frequency.setText(str(self.center.frequency))
	self.controls.entry_spindown.setText(str(self.center.spindown))
	self.controls.entry_ra.setText(str(self.center.ra))
	self.controls.entry_dec.setText(str(self.center.dec))
	self.controls.entry_iota.setText(str(self.center.iota))
	self.controls.entry_psi.setText(str(self.center.psi))

	self.controls.entry_step_frequency.setText(str(self.center.frequency_step))
	self.controls.entry_step_spindown.setText(str(self.center.spindown_step))
	self.controls.entry_step_ra.setText(str(self.center.ra_step))
	self.controls.entry_step_dec.setText(str(self.center.dec_step))
	self.controls.entry_step_iota.setText(str(self.center.iota_step))
	self.controls.entry_step_psi.setText(str(self.center.psi_step))

    def get_controls(self):
	self.center.frequency=float(self.controls.entry_frequency.text())
	self.center.spindown=float(self.controls.entry_spindown.text())
	self.center.ra=float(self.controls.entry_ra.text())
	self.center.dec=float(self.controls.entry_dec.text())
	self.center.iota=float(self.controls.entry_iota.text())
	self.center.psi=float(self.controls.entry_psi.text())

	self.center.frequency_step=float(self.controls.entry_step_frequency.text())
	self.center.spindown_step=float(self.controls.entry_step_spindown.text())
	self.center.ra_step=float(self.controls.entry_step_ra.text())
	self.center.dec_step=float(self.controls.entry_step_dec.text())
	self.center.iota_step=float(self.controls.entry_step_iota.text())
	self.center.psi_step=float(self.controls.entry_step_psi.text())

    def initMainWidget(self):
	#self.layout1 = QGridLayout()
	
	self.f=QWidget(self)
	self.setCentralWidget(self.f)
	
	self.layout1=QGridLayout(self.f)
	
	self.controls = QWidget(self)
	self.layout1.addWidget(self.controls, 1, 2)
	
	self.controls.layout2 = QGridLayout(self.controls)
	
	self.controls.display=QLabel("", self.controls)
	self.controls.layout2.addWidget(self.controls.display, 0, 0, 1, 3)
	
	
	self.controls.header_var=QLabel("", self.controls)
	self.controls.header_value=QLabel("Value", self.controls)
	self.controls.header_step=QLabel("Step", self.controls)
	
	self.controls.layout2.addWidget(self.controls.header_var, 1, 0)
	self.controls.layout2.addWidget(self.controls.header_value, 1, 1)
	self.controls.layout2.addWidget(self.controls.header_step, 1, 2)
	
	self.controls.label_frequency = QLabel("Frequency", self.controls)
	self.controls.label_spindown = QLabel("Spindown", self.controls)
	self.controls.label_ra = QLabel("RA", self.controls)
	self.controls.label_dec = QLabel("DEC", self.controls)
	self.controls.label_iota = QLabel("iota", self.controls)
	self.controls.label_psi = QLabel("psi", self.controls)

	self.controls.layout2.addWidget(self.controls.label_frequency, 2, 0)
	self.controls.layout2.addWidget(self.controls.label_spindown, 3, 0)
	self.controls.layout2.addWidget(self.controls.label_ra, 4, 0)
	self.controls.layout2.addWidget(self.controls.label_dec, 5, 0)
	self.controls.layout2.addWidget(self.controls.label_iota, 6, 0)
	self.controls.layout2.addWidget(self.controls.label_psi, 7, 0)

	self.controls.entry_frequency = QLineEdit("Frequency", self.controls)
	self.controls.entry_spindown = QLineEdit("Spindown", self.controls)
	self.controls.entry_ra = QLineEdit("RA", self.controls)
	self.controls.entry_dec = QLineEdit("DEC", self.controls)
	self.controls.entry_iota = QLineEdit("iota", self.controls)
	self.controls.entry_psi = QLineEdit("psi", self.controls)

	self.controls.layout2.addWidget(self.controls.entry_frequency, 2, 1)
	self.controls.layout2.addWidget(self.controls.entry_spindown, 3, 1)
	self.controls.layout2.addWidget(self.controls.entry_ra, 4, 1)
	self.controls.layout2.addWidget(self.controls.entry_dec, 5, 1)
	self.controls.layout2.addWidget(self.controls.entry_iota, 6, 1)
	self.controls.layout2.addWidget(self.controls.entry_psi, 7, 1)

	self.controls.entry_step_frequency = QLineEdit("Frequency", self.controls)
	self.controls.entry_step_spindown = QLineEdit("Spindown", self.controls)
	self.controls.entry_step_ra = QLineEdit("RA", self.controls)
	self.controls.entry_step_dec = QLineEdit("DEC", self.controls)
	self.controls.entry_step_iota = QLineEdit("iota", self.controls)
	self.controls.entry_step_psi = QLineEdit("psi", self.controls)

	self.controls.layout2.addWidget(self.controls.entry_step_frequency, 2, 2)
	self.controls.layout2.addWidget(self.controls.entry_step_spindown, 3, 2)
	self.controls.layout2.addWidget(self.controls.entry_step_ra, 4, 2)
	self.controls.layout2.addWidget(self.controls.entry_step_dec, 5, 2)
	self.controls.layout2.addWidget(self.controls.entry_step_iota, 6, 2)
	self.controls.layout2.addWidget(self.controls.entry_step_psi, 7, 2)
	
	self.controls.entry_frequency.setMinimumWidth(100)
	self.controls.entry_step_frequency.setMinimumWidth(100)

	self.controls.recompute = QPushButton("Recompute", self.controls)
	self.controls.layout2.addWidget(self.controls.recompute, 8, 0, 1, 3)

	self.controls.best_candidate = QPushButton("Best candidate", self.controls)
	self.controls.layout2.addWidget(self.controls.best_candidate, 9, 0, 1, 3)

	self.set_controls()
	
	self.disp00=two_display(self, "frequency", "spindown")
	self.layout1.addWidget(self.disp00, 0, 0)
	
	self.disp01=two_display(self, "ra", "dec")
	self.layout1.addWidget(self.disp01, 0, 1)

	self.disp10=two_display(self, "spindown", "dec")
	self.layout1.addWidget(self.disp10, 1, 0)

	self.disp11=two_display(self, "spindown", "ra")
	self.layout1.addWidget(self.disp11, 1, 1)

	self.disp02=two_display(self, "iota", "psi")
	self.layout1.addWidget(self.disp02, 0, 2)

	#self.update(self.center)

    def best_candidate(self):
	    c=None
	    for attr in [ "00", "01", "02", "10", "11" ] :
		    w=getattr(self, "disp"+attr)
		    c2=w.best_candidate()
		    if c2 == None : continue
		    if c==None :
			    c=c2
			    continue
		    if c2.snr>c.snr : c=c2
	    if c==None : return
	    self.center.copy_coords(c)
	    c=compute(self.center)
	    self.controls.display.setText("snr: %f\nstrain: %g\npower_cor: %f\nifo_freq: %f\n" % (c.snr, c.strain, c.power_cor, c.ifo_freq))
	    self.set_controls()

    def recompute(self):
	    self.get_controls()
	    self.update(self.center)

    def update(self, center):
	    self.disp00.update(center)
	    self.disp01.update(center)
	    self.disp02.update(center)
	    self.disp10.update(center)
	    self.disp11.update(center)
	    #self.disp12.update(center)
	

def main(args):
    app=QApplication(args)
    mainWindow = powerflux(app)
    mainWindow.show()
    app.connect(app, SIGNAL("lastWindowClosed()"), app, SLOT("quit()"))
    app.connect(mainWindow.controls.recompute, SIGNAL("clicked()"), mainWindow.recompute)
    app.connect(mainWindow.controls.best_candidate, SIGNAL("clicked()"), mainWindow.best_candidate)
    #mainWindow.update(mainWindow.center)
    app.exec_()



if __name__ == "__main__":
    main(sys.argv)

