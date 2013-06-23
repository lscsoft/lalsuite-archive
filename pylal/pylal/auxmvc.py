# Module to keep definitions of job and node classes for auxmvc pipeline. 

import os
import sys
import tempfile
from glue import pipeline
import tempfile

def construct_command(node):
  command_string = node.job().get_executable() + " " + node.get_cmd_line()
  return [opt for opt in command_string.split(" ") if opt.strip()]
  

class auxmvc_DAG(pipeline.CondorDAG):
        def __init__(self, basename, log_path):
                self.basename = basename
                tempfile.tempdir = log_path
                tempfile.template = self.basename + '.dag.log.'
                logfile = tempfile.mktemp()
                fh = open(logfile, "w" )
                fh.close()
                pipeline.CondorDAG.__init__(self,logfile)
                self.set_dag_file(self.basename)
                self.jobsDict = {}
                #self.id = 0
        def add_node(self, node):
                #self.id+=1
                pipeline.CondorDAG.add_node(self, node)


  
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
    print exec_name
    print executable
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
 
  
  
class build_auxmvc_vectors_job(auxmvc_analysis_job):
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
  
class prepare_training_auxmvc_samples_job(auxmvc_analysis_job):
  """
  Job for preparing training auxmvc samples. 
  """
  def __init__(self, cp):
	"""
	"""
	sections = ['prepare-training-auxmvc-samples']
	exec_name = 'idq_prepare_training_auxmvc_samples'
	tag_base = 'training_auxmvc'
	auxmvc_analysis_job.__init__(self,cp,sections,exec_name,tag_base=tag_base)	
		
class prepare_training_auxmvc_samples_node(pipeline.CondorDAGNode):
  """
  Node for preparing training auxmvc samples job. 
  """
  def __init__(self, job, source_dir, gps_start_time, gps_end_time, output_file, dq_segments="", dq_segments_name="",p_node=[]):
    job.set_stdout_file('logs/' + output_file.replace('.pat', '.out'))
    job.set_stderr_file('logs/' + output_file.replace('.pat', '.err'))
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_output_file(output_file)
    self.add_opt('source-directory', source_dir)
    self.add_opt('gps-start-time', gps_start_time)
    self.add_opt('gps-end-time', gps_end_time)
    if dq_segments and dq_segments_name:
	self.add_opt('dq-segments', dq_segments)
	self.add_opt('dq-segments-name', dq_segments_name)
    self.add_opt('output-file', output_file)
    for p in p_node:
      self.add_parent(p)

    

class train_forest_job(auxmvc_analysis_job):
  """
  Training job for random forest (MVSC). 
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
  def __init__(self, job, training_data_file, p_node=[]):
    pipeline.CondorDAGNode.__init__(self,job)
    self.job.set_stdout_file('logs/' + training_data_file.replace('.pat', '.out'))
    self.job.set_stderr_file('logs/' + training_data_file.replace('.pat', '.err'))
    self.add_input_file(training_data_file)
    self.training_data_file = self.get_input_files()[0]
    self.trainedforest = self.training_data_file.replace('.pat','.spr')
	self.add_output_file(self.trainedforest)	
	self.add_short_opt("f", self.trainedforest)
    self.add_file_arg(" %s" % (self.training_data_file))
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
		

########################  svm for idq  ############################


class create_DAG(pipeline.CondorDAG):
    """
    """
    def __init__(self, config_file, log_path, run_tag):
        self.config_file = str(config_file)
        self.basename = self.config_file.replace('.ini','')+run_tag
        tempfile.tempdir = log_path
        tempfile.template = self.basename + '.dag.log.'
        logfile = tempfile.mktemp()
        fh = open( logfile, "w" )
        fh.close()
        pipeline.CondorDAG.__init__(self,logfile)
	self.set_dag_file(self.basename)
        self.jobsDict = {}
        self.id = 0

    def add_node(self, node):
        self.id+=1
        pipeline.CondorDAG.add_node(self, node)


class use_svm_job(auxmvc_analysis_job):
    """
    """
    def __init__(self, cp):
        """
        """
        sections = []
        exec_name = 'svm_evaluate_cmd'
        tag_base  = 'svm_evaluate'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name, tag_base=tag_base, short_opts=True)


class use_svm_node(pipeline.CondorDAGNode):
    """
    Node for SVM evaluation job.
    """
    def __init__(self, job, cp, test_file, range_file, svm_model, predict_file, p_node=[]):
        job.set_stdout_file('logs/' + predict_file.replace('.dat', '.out'))
        job.set_stderr_file('logs/' + predict_file.replace('.dat', '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(test_file)
        self.add_input_file(range_file)
        self.add_input_file(svm_model)
        self.add_output_file(predict_file)

        self.scale_cmd = cp.get('svm_evaluate','svm_scale_cmd')
        self.predict_cmd = cp.get('svm_evaluate', 'svm_predict_cmd')
	self.test_file = self.get_input_files()[0]
        self.range_file = self.get_input_files()[1]
        self.svm_model = self.get_input_files()[2]
        self.predict_file = self.get_output_files()[0]
	self.add_file_arg(" --scale %s --predict %s -i %s -r %s -m %s -o %s" % (self.scale_cmd, self.predict_cmd, self.test_file, self.range_file, self.svm_model, self.predict_file))
        for p in p_node:
            self.add_parent(p)


class train_svm_job(auxmvc_analysis_job):
  """
  Training job for svm.
  """
  def __init__(self, cp):
    """
    """
    sections = [] #not section in configuration yet
    exec_name = 'svm_train_cmd'
    tag_base  = 'svm_train'
    auxmvc_analysis_job.__init__(self, cp, sections, exec_name, tag_base=tag_base, short_opts=True)


class train_svm_node(pipeline.CondorDAGNode):
    """
    Node for SVM train job.
    """
    def __init__(self, job, dag, cp, train_file, range_file, model_file, p_node=[]):
        job.set_stdout_file('logs/' + train_file.replace('.dat', '.out'))
        job.set_stderr_file('logs/' + train_file.replace('.dat', '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(train_file)
        self.add_output_file(range_file)
        self.add_output_file(model_file)
        self.scale_cmd = cp.get('svm_evaluate','svm_scale_cmd')
        self.train_cmd = cp.get('svm_evaluate','svm_train_cmd')
        self.gamma = cp.get('svm_evaluate','svm_gamma')
        self.cost  = cp.get('svm_evaluate','svm_cost')
        self.train_file = self.get_input_files()[0]
        self.range_file = self.get_output_files()[0]
        self.model_file = self.get_output_files()[1]

        self.train_file_svm = os.path.abspath(self.train_file) + '.mid'
        self.scale_file = os.path.abspath(self.train_file) + '.scale'
        self.add_file_arg(" --scale %s --train %s --train-file %s --train-file-svm %s --scale-file %s --range-file %s --model-file %s -g %s -c %s " % (self.scale_cmd, self.train_cmd, self.train_file, self.train_file_svm, self.scale_file, self.range_file, self.model_file, self.gamma, self.cost))
        self.set_post_script("/bin/rm ")
        self.add_post_script_arg(self.train_file_svm)
        self.add_post_script_arg(self.scale_file)

        for p in p_node:
            self.add_parent(p)
	dag.add_node(self)

