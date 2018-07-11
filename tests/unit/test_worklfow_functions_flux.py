#################################################################################################################
@pytest.mark.skip(reason="only on FLUX")
def test_submit_flux_job_one_job():
    output_directory = "/scratch/hmobley_fluxod/annasint/code/tests/test_data/"
    suffix = "test_script1"
    today = dt.datetime.today().strftime("%Y-%m-%d")  # todo refactor today variable
    job_name = "test_job"
    script = "mkdir TEST1"
    output = workflow.submit_flux_job(output_directory, suffix, today,
                                      job_name, script)
    assert type(int(output)) == int


@pytest.mark.skip(reason="only on FLUX")
def test_submit_flux_job_two_jobs():
    # Submit first job:
    output_directory = "/scratch/hmobley_fluxod/annasint/code/tests/test_data/"
    suffix = "test_script1"
    today = dt.datetime.today().strftime("%Y-%m-%d")
    job_name = "test_job"
    script = "mkdir TEST1"

    output1 = workflow.submit_flux_job(output_directory, suffix, today,
                                       job_name, script)

    # Submit second job:
    suffix2 = 'test_script2'
    script2 = "cd TEST1\nmkdir TEST2"

    output2 = workflow.submit_flux_job(output_directory, suffix2, today,
                                       job_name, script2, job_dependency=output1)
    assert type(int(output2)) == int

# todo make tests clean up after themselves

# this will test both on flux, plus if job dependency works properly


@pytest.mark.skip(reason="only on FLUX")
def test_run_fastqc_after_run_trim_job():
    fastq_file_input = "/scratch/hmobley_fluxod/annasint/code/data/reads/HM86_UR_copy.fastq.gz"
    output_directory = "/scratch/hmobley_fluxod/annasint/code/tests/test_data"
    today = dt.datetime.today().strftime("%Y-%m-%d")
    config_dict = workflow.process_config(config_file="config")
    local = False
    output_file_name, jobid = workflow.run_trim_job(fastq_file_input, output_directory,
                                                    today, config_dict, local)
    assert type(int(jobid)) == int  # only way can figure out if the job has been submitted
    actual_out_dir, jobid2 = workflow.run_fastqc_job(output_file_name, output_directory, today,
                                                     config_dict, local, job_dependency=jobid)
    assert type(int(jobid2)) == int  # todo come up with a better assert statement



