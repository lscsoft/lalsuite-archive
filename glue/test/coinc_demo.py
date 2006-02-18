import math
import random

from glue.ligolw import ligolw
from glue.ligolw import lsctables

processids = iter(lsctables.ProcessIDs())
burstids = iter(lsctables.SnglBurstIDs())
inspids = iter(lsctables.SnglInspiralIDs())
slideids = iter(lsctables.TimeSlideIDs())
coincids = iter(lsctables.CoincIDs())

def new_proc(prog):
	row = lsctables.Process()
	row.process_id = str(processids.next())
	row.program = prog
	return row

def new_burst(process_id):
	row = lsctables.SnglBurst()
	row.process_id = process_id
	row.peak_time = int(random.random() * 900000000)
	row.peak_time_ns = int(random.random() * 1000000000)
	row.confidence = random.random() * -1000.0
	row.event_id = str(burstids.next())
	return row

def new_insp(process_id):
	row = lsctables.SnglInspiral()
	row.process_id = process_id
	row.end_time = int(random.random() * 900000000)
	row.end_time_ns = int(random.random() * 1000000000)
	row.snr = random.random() * 100.0
	row.event_id = str(inspids.next())
	return row

def new_slide(process_id, h1off, h2off, l1off):
	h1 = lsctables.TimeSlide()
	h2 = lsctables.TimeSlide()
	l1 = lsctables.TimeSlide()

	h1.instrument = "H1"
	h2.instrument = "H2"
	l1.instrument = "L1"

	h1.time_slide_id = h2.time_slide_id = l1.time_slide_id = str(slideids.next())
	h1.process_id = h2.process_id = l1.process_id = process_id

	h1.offset = h1off
	h2.offset = h2off
	l1.offset = l1off

	return h1, h2, l1

def new_coinc(process_id, time_slide_id):
	row = lsctables.Coinc()
	row.process_id = process_id
	row.coinc_event_id = str(coincids.next())
	row.time_slide_id = time_slide_id
	return row

def new_coinc_map(coinc_event_id, event_id):
	row = lsctables.CoincMap()
	row.coinc_event_id = coinc_event_id
	row.event_id = event_id
	return row


doc = ligolw.Document()
llw = ligolw.LIGO_LW()
doc.appendChild(llw)

process = lsctables.New(lsctables.ProcessTable, ["process_id", "program"])
llw.appendChild(process)
burst_proc = new_proc("lalapps_power")
insp_proc = new_proc("lalapps_inspiral")
slide_proc = new_proc("lalapps_pick_slides")
coinc_proc = new_proc("lalapps_ubercoinc")
process.append(burst_proc)
process.append(insp_proc)
process.append(slide_proc)
process.append(coinc_proc)

burst = lsctables.New(lsctables.SnglBurstTable, ["process_id", "peak_time", "peak_time_ns", "confidence", "event_id"])
llw.appendChild(burst)
burst.append(new_burst(burst_proc.process_id))
burst.append(new_burst(burst_proc.process_id))
burst.append(new_burst(burst_proc.process_id))

insp = lsctables.New(lsctables.SnglInspiralTable, ["process_id", "end_time", "end_time_ns", "snr", "event_id"])
llw.appendChild(insp)
insp.append(new_insp(insp_proc.process_id))
insp.append(new_insp(insp_proc.process_id))
insp.append(new_insp(insp_proc.process_id))

slide = lsctables.New(lsctables.TimeSlideTable)
llw.appendChild(slide)
map(slide.append, new_slide(slide_proc.process_id, 0.0, 0.0, 0.0))
map(slide.append, new_slide(slide_proc.process_id, 0.0, 0.0, 50.0))

coinc = lsctables.New(lsctables.CoincTable)
llw.appendChild(coinc)
coinc.append(new_coinc(coinc_proc.process_id, slide[0].time_slide_id))
coinc.append(new_coinc(coinc_proc.process_id, slide[0].time_slide_id))

coincmap = lsctables.New(lsctables.CoincMapTable)
llw.appendChild(coincmap)
coincmap.append(new_coinc_map(coinc[0].coinc_event_id, burst[0].event_id))
coincmap.append(new_coinc_map(coinc[0].coinc_event_id, burst[1].event_id))
coincmap.append(new_coinc_map(coinc[0].coinc_event_id, insp[0].event_id))
coincmap.append(new_coinc_map(coinc[1].coinc_event_id, burst[2].event_id))
coincmap.append(new_coinc_map(coinc[1].coinc_event_id, insp[1].event_id))

doc.write()
