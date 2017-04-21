from seq_jobs import do_jobs

do_jobs(command = "sh run.sh",
        files = ["one_sto.json", "run.sh"],
        out_dir = "out",
        label_vals = [("DT",     [0.1, 0.01]),
                      ("METHOD", ["Euler", "RK4"])],
        dir_name_rule = "default",
        log_file = "seq_jobs.log"
)
