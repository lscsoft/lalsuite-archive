# Module to keep definitions of job and node classes for auxmvc pipeline. 
import sys
from glue import pipeline

def construct_command(node):
  command_string = node.job().get_executable() + " " + node.get_cmd_line()
  return [opt for opt in command_string.split(" ") if opt.strip()]
  
  
def build_auxmvc_vectors(triggger_dict, main_channel, time_window, signif_threshold, clean_samples_rate, output_file_name, channels=None, unsafe_channels=None):
  """
  Given dictionary of triggers from multiple channels, the function constructs auxmvc vectors out of them. Result is saved into output_file_name. 
  """  

  # use only channels from the channels file, if provided
  if channels:
    selected_channels = event.read_channels_from_file(channels)
    trigger_dict.keep_channels(selected_channels)
  
    # to ensure consistency in dimensionality of auxmvc vectors
    # add the channels from the selected channels that are absent in trigger_dict
    for channel in selected_channels: 
	  if not channel in trigger_dict.channels():
	    trigger_dict[channel] = []


  # get rid of unsafe channels if provided
  if unsafe_channels:
    unsafe_channels = event.read_channels_from_file(unsafe_channels)
    trigger_dict.remove_channels(unsafe_channels)
  

  # keep only the triggers from the [gps_start_time, gps_end_time] segment
  if opts.gps_start_time and opts.gps_end_time:
  
    #first keep all triggers from the segment expanded by the time concidence window, so that not to loose coincidences
	trigger_dict.include([[gps_start_time - time_window, gps_end_time + time_window]])
	#then in the main channel keep triggers that fall within the segment.
	trigger_dict.include([[gps_start_time, gps_end_time]], channels=[main_channel])

  # apply significance threshold to the triggers from the main channel
  trigger_dict.apply_signif_threshold(channels = [main_channel], threshold = signif_threshold)

  # construct glitch auxmvc vectors
  aux_glitch_vecs = event.build_auxmvc_vectors(trigger_dict, main_channel = main_channel, coincidence_time_window = time_window)

  # compute number of clean samples
  n_clean = round(clean_samples_rate * (gps_end_time - gps_start_time))

  # generate random times for clean samples
  clean_times = numpy.random.uniform(gps_start_time, gps_end_time, n_clean)

  # construct clean auxmvc vectors
  aux_clean_vecs = event.build_auxmvc_vectors(trigger_dict, main_channel = main_channel, coincidence_time_window = time_window, build_clean_samples = True, clean_times = clean_times)

  # get rid of clean samples that are near real triggers in the main channel.
  aux_clean_vecs = auxmvc_utils.get_clean_samples(aux_clean_vecs)


  # convert glitch and clean auxmvc vectors into MVSC evaluation set
  mvsc_evaluation_set=auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers = aux_glitch_vecs, KWAuxCleanTriggers = aux_clean_vecs)
	  
  # save MVSC evaluation set in file	  

  auxmvc_utils.WriteMVSCTriggers(mvsc_evaluation_set, output_filename = output_file_name, Classified = False)
  
  return mvsc_evaluation_set

  
  
  
  
#####################  JOB and NODE classes for auxmvc pipeline  #################################  
  
  
class auxmvc_analysis_job(pipeline.AnalysisJob, pipeline.CondorDAGJob):
  """
  A basic auxmvc job class. Sets common atributes needed for any auxmvc job. It uses config parser object to 
  set the options. 
  """
  def __init__(self,cp,sections,exec_name,tag_base='', id ='',extension='',dax=False, short_opts=False):
    """
    cp = ConfigParser object from which options are read.
    sections = sections of the ConfigParser that get added to the opts
    exec_name = exec_name name in ConfigParser
    """
    self.__exec_name = exec_name
    self.__extension = extension
    self.tag_base = tag_base
    universe = cp.get('condor','universe')
    executable = cp.get('condor',exec_name)
    pipeline.CondorDAGJob.__init__(self,universe,executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)
    self.add_condor_cmd('copy_to_spool','False')
    self.add_condor_cmd('getenv','True')
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.__use_gpus = cp.has_option('condor', 'use-gpus')
    for sec in sections:
      if cp.has_section(sec):
        if short_opts:
          self.add_short_ini_opts(cp,sec)
        else:
          self.add_ini_opts(cp, sec)
      else:
        print >>sys.stderr, "warning: config file is missing section [" + sec + "]"

    self.set_stdout_file('logs/' + tag_base + id + '.out')
    self.set_stderr_file('logs/' + tag_base + id + '.err')
    self.set_sub_file(tag_base + '.sub')

  def set_exec_name(self,exec_name):
    """
    Set the exec_name name
    """
    self.__exec_name = exec_name
	
  def set_exec_name(self,exec_name):
    """
    Set the exec_name name
    """
    self.__exec_name = exec_name

  def get_exec_name(self):
    """
    Get the exec_name name
    """
    return self.__exec_name

  def set_extension(self,extension):
    """
    Set the file extension
    """
    self.__extension = extension

  def get_extension(self):
    """
    Get the extension for the file name
    """
    return self.__extension

  def get_use_gpus(self):
    """
    Get whether this job was requested to run on a GPU node
    """
    return self.__use_gpus
	
  def add_short_ini_opts(self, cp, section):
    """
    Parse command line options from a given section in an ini file and
    pass to the executable as short options.
    @param cp: ConfigParser object pointing to the ini file.
    @param section: section of the ini file to add to the options.
    """
    for opt in cp.options(section):
      arg = cp.get(section,opt)
      self.add_short_opt(opt,arg)
 
  
  
class build_auxmvc_vectors_job(pipeline.CondorDAGJob):
  """
  Job for building auxmvc feature vectors. 
  """
  def __init__(self, cp):
	"""
	"""
	sections = ['build-auxmvc-vectors']
	exec_name = 'idq_build_auxmvc_vectors'
	tag_base = 'build_vectors'
	auxmvc_analysis_job.__init__(self,cp,sections,exec_name,tag_base=tag_base)	
		

class build_auxmvc_vectors_node(pipeline.CondorDAGNode):
  """
  Dag node for building axumvc feature vector.
  """
  def __init__(self, job, triggerfiles, gps_start_time, gps_end_time, p_node=[]):
    pipeline.CondorDAGNode.__init__(self,job)
    self.job.set_stdout_file('logs/' + job.tag_base + '-' + str(gps_start_time) + '-' + str(gps_end_time - gps_start_time) +'.out')
    self.job.set_stderr_file('logs/' + job.tag_base + '-' +  str(gps_start_time) + '-' + str(gps_end_time - gps_start_time) +'.err')
    self.add_input_file(triggerfiles)
    trigger_type = self.job.get_opt('trigger-type')
    ifo = self.job.get_opt('ifo')
    user_tag = self.job.get_opt('user-tag')
    output_file_name= trigger_type + ifo + "-" + "evaluation_" + user_tag  + "-" + gps_start_time + "-" + str(gps_end_time - gps_start_time) + ".pat"
    self.add_output_file(output_file_name)
    self.add_var_opt('trigger-files', triggerfiles)
    for p in p_node:
      self.add_parent(p)
  

class train_forest_job(pipeline.CondorDAGJob):
  """
  Training job for random forets (MVSC). 
  """
  def __init__(self, cp):
    """
    """
    sections = ['train_forest']
    exec_name = 'SprBaggerDecisionTreeApp'
    tag_base = 'train_forest'
    auxmvc_analysis_job.__init__(self,cp,sections,exec_name,tag_base=tag_base, short_opts=True)	
		

class train_forest_node(pipeline.CondorDAGNode):
  """
  Dag node for training the random forest (MVSC).
  """
  def __init__(self, job, trainingdatafile, p_node=[]):
    pipeline.CondorDAGNode.__init__(self,job)
    self.job.set_stdout_file('logs/' + trainingdatafile.replace('.pat', '.out'))
    self.job.set_stderr_file('logs/' + trainingdatafile.replace('.pat', '.err'))
    self.add_input_file(trainingdatafile)
    self.trainingdatafile = self.get_input_files()[0]
    self.trainedforest = self.trainingdatafile.replace('.pat','.spr')
    #self.trainedout = self.trainingdatafile.replace('.pat','.out')
    self.add_file_arg(" %s %s" % (self.trainedforest, self.trainingdatafile))
    self.add_output_file(self.trainedforest)
    for p in p_node:
      self.add_parent(p)

class use_forest_job(auxmvc_analysis_job):
  """
  Job using random forest to evaluate unclassified data.
  """
  def __init__(self, cp):
    """
    """
    sections = ['forest_evaluate']
    exec_name = 'SprOutputWriterApp'
    tag_base = 'forest_evaluate'
    auxmvc_analysis_job.__init__(self,cp, sections, exec_name, tag_base=tag_base, short_opts=True)
    self.add_short_opt("A", "")

class use_forest_node(pipeline.CondorDAGNode):
  """
  Node for radnom forest evaluation job. 
  """
  def __init__(self, job, trainedforest, file_to_rank, ranked_file,p_node=[]):
    job.set_stdout_file('logs/' + ranked_file.replace('.dat', '.out'))
    job.set_stderr_file('logs/' + ranked_file.replace('.dat', '.err'))
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_input_file(trainedforest)
    self.add_input_file(file_to_rank)
    self.add_output_file(ranked_file)
    self.trainedforest = self.get_input_files()[0]
    self.file_to_rank = self.get_input_files()[1]
    self.ranked_file = ranked_file
    self.add_file_arg(" %s %s %s" % (self.trainedforest, self.file_to_rank, self.ranked_file))
    for p in p_node:
      self.add_parent(p)


class forest_add_excluded_vars_job(auxmvc_analysis_job):
  """
  A simple fix job that adds the variables excluded by forest (MVSC) from classification into the output file.
  Need to be run right after use_forest job. 
  """
  def __init__(self, cp):
    """
    """
    sections = ['forest_add_excluded_vars']
    exec_name = 'forest_add_excluded_vars'
    tag_base = 'forest_add_excluded_vars'
    auxmvc_analysis_job.__init__(self,cp, sections, exec_name, tag_base=tag_base)
    self.add_opt('excluded-variables', cp.get('forest_evaluate', 'z'))
class forest_add_excluded_vars_node(pipeline.CondorDAGNode):
  """
  Node for forest_add_excluded_vars_job.
  """
  def __init__(self, job, patfile, datfile, p_node=[]):
    job.set_stdout_file('logs/' + datfile.replace('.dat', 'faev.out'))
    job.set_stderr_file('logs/' + datfile.replace('.dat', 'faev.err'))
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_input_file(patfile)
    self.add_input_file(datfile)
    self.add_var_opt('pat-file',patfile)
    self.add_var_opt('dat-file',datfile)
    for p in p_node:
      self.add_parent(p)




class plot_channels_significance_job(pipeline.CondorDAGJob):
  """   
  Job that makes verious plots and histograms using significance of the axuiloary channels.      
  """
  def __init__(self, cp):
    sections = ['plot-forest-channels-significance']
    exec_name = 'auxmvc_plot_mvsc_channels_significance'
    tag_base = 'plot_channels_signif'
    auxmvc_analysis_job.__init__(self,cp, sections, exec_name, tag_base=tag_base)      
				
	
class plot_channels_significance_node(pipeline.CondorDAGNode):
  """
  Node for plot_channels_significance_job.
  """
  def __init__(self, job, input, p_node=[]):
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_opt("input", input)
    for p in p_node:                       
      self.add_parent(p)
			  


class result_plots_job(pipeline.CondorDAGJob):
  """
  Job that makes plots based on results of evaluation e.g. ROC curves.
  """
  def __init__(self, cp, tag_base='RESULT_PLOTS'):
    """
    """
    sections = ['result_plots']
    exec_name = 'auxmvc_result_plots'
    tag_base = 'auxmvc_result_plots'
    auxmvc_analysis_job.__init__(self,cp, sections, exec_name, tag_base=tag_base)  


class result_plots_node(pipeline.CondorDAGNode):
  """
  Node for result_plots_job.
  """
  def __init__(self, job, datfiles, p_node=[]):
    pipeline.CondorDAGNode.__init__(self,job)
    for file in datfiles:
      self.add_file_arg(file[0])
    for p in p_node:
      self.add_parent(p)
		


