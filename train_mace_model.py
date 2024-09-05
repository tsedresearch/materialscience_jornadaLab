import subprocess

command = [
    "python3", "/global/homes/t/tsedalya/env_folder/allegro_env/lib/python3.9/site-packages/mace/cli/run_train.py",
    "--name", "MACE_model",
    "--train_file", "/{the training file directory}",
    "--test_file", "/{the testing file directory}",
    "--config_type_weights", '{"Default":1.0}',
    "--E0s", 'average',
    "--model", "MACE",
    "--hidden_irreps", '128x0e + 128x1o',
    "--r_max", "5.0",
    "--batch_size", "10",
    "--max_num_epochs", "1000",
    "--swa",
    "--start_swa", "1200",
    "--ema",
    "--ema_decay", "0.99",
    "--amsgrad",
    "--device", "cuda",
    "--model_dir", "/{the directory to store the trained model}",
    "--results_dir", "/{the directory to store the results of the training of the model}",
    "--log_dir", "/{the directory to store the log file}
     ]

subprocess.run(command)
