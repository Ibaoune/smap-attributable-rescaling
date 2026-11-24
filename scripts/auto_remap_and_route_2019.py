#!/usr/bin/env python3
"""
Automated Orchestrator for 2019 DA experiments:
1. Monitors DA_nocdf_2019, DA_cdftrain_2019, and DA_v2_2019 chains.
2. When all three reach 2020-01-01 without errors:
   - Runs the three runoff remapping scripts.
   - Submits the three HyMAP2 2019 routing jobs on Toubkal (compute partition, empowermed-ahl6xm8o7mg-DEFAULT-CPU).
   - Monitors HyMAP runs until completion.
   - Checks and verifies the existence of streamflow files in routing/output/P3_2019/<EXP>/ROUTING/.
   - Generates status report (without calculating any Decision 16 criteria).
"""

import os
import sys
import time
import subprocess
import glob
from datetime import datetime

WORKSPACE = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling"
P1_DIR = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco"

EXPERIMENTS = [
    "DA_nocdf_2019",
    "DA_cdftrain_2019",
    "DA_v2_2019",
]

def log(msg):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {msg}", flush=True)

def check_chain_completed(exp):
    rst_file = os.path.join(WORKSPACE, "experiments/da_v2", exp, "restarts/surf/LIS_RST_NOAHMP401_202001010000.d01.nc")
    return os.path.isfile(rst_file)

def get_active_jobs():
    cmd = ["squeue", "-u", "mohammad.elaabaribao", "-h", "-o", "%i %j %T"]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    jobs = {}
    for line in res.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 3:
            jobs[parts[0]] = {"name": parts[1], "state": parts[2]}
    return jobs

def check_for_errors(exp):
    log_dir = os.path.join(WORKSPACE, "experiments/da_v2", exp, "logs")
    err_files = glob.glob(f"{log_dir}/**/*.err", recursive=True)
    for ef in err_files:
        if os.path.getsize(ef) > 0:
            with open(ef, "r") as f:
                content = f.read()
            if "MPI_ABORT" in content or "ERROR" in content or "Error opening" in content:
                return ef, content
    return None, None

def run_remapping():
    log("=== Commencing Runoff Remapping for all 3 experiments ===")
    scripts = {
        "DA_nocdf_2019": os.path.join(P1_DIR, "experiments/NorthMor/routing/scripts/remap_runoff_da_nocdf_2019.py"),
        "DA_cdftrain_2019": os.path.join(P1_DIR, "experiments/NorthMor/routing/scripts/remap_runoff_da_cdftrain_2019.py"),
        "DA_v2_2019": os.path.join(P1_DIR, "experiments/NorthMor/routing/scripts/remap_runoff_da_v2_2019.py"),
    }
    
    for exp, script_path in scripts.items():
        log(f"Running remapping for {exp} using {script_path}...")
        t0 = time.time()
        cmd = ["python3", script_path]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        elapsed = time.time() - t0
        if res.returncode != 0:
            log(f"ERROR: Remapping failed for {exp} with returncode {res.returncode}:\n{res.stderr}")
            raise RuntimeError(f"Remapping failed for {exp}")
        log(f"Remapping completed for {exp} in {elapsed:.1f}s.")
        
        # Verify output files exist
        remap_dir = os.path.join(P1_DIR, "experiments/NorthMor/routing/runoff_remap_180x120", exp, "SURFACEMODEL")
        files = glob.glob(f"{remap_dir}/**/*.nc", recursive=True)
        log(f"Remapped files count for {exp}: {len(files)} files found.")
        if len(files) < 365:
            log(f"WARNING: Expected at least 365 files, found {len(files)}.")

def submit_hymap_jobs():
    log("=== Submitting HyMAP2 2019 Routing Jobs ===")
    jobs = {
        "DA_nocdf_2019": os.path.join(P1_DIR, "experiments/NorthMor/routing/jobs/run_hymap2_da_nocdf_2019.sh"),
        "DA_cdftrain_2019": os.path.join(P1_DIR, "experiments/NorthMor/routing/jobs/run_hymap2_da_cdftrain_2019.sh"),
        "DA_v2_2019": os.path.join(P1_DIR, "experiments/NorthMor/routing/jobs/run_hymap2_da_v2_2019.sh"),
    }
    
    submitted_ids = {}
    for exp, job_script in jobs.items():
        cmd = ["sbatch", job_script]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            log(f"ERROR submitting HyMAP job for {exp}: {res.stderr}")
            raise RuntimeError(f"Failed to submit HyMAP job for {exp}")
        
        # Parse job id: "Submitted batch job 123456"
        job_id = res.stdout.strip().split()[-1]
        submitted_ids[exp] = job_id
        log(f"Submitted HyMAP job for {exp}: Job ID {job_id}")
    
    return submitted_ids

def monitor_hymap_jobs(submitted_ids):
    log("=== Monitoring HyMAP2 Jobs ===")
    remaining = dict(submitted_ids)
    while remaining:
        active = get_active_jobs()
        for exp, jid in list(remaining.items()):
            if jid not in active:
                log(f"HyMAP job for {exp} (Job ID {jid}) has finished.")
                del remaining[exp]
            else:
                log(f"HyMAP job for {exp} (Job ID {jid}) is {active[jid]['state']}...")
        if remaining:
            time.sleep(30)
    log("All HyMAP jobs have completed.")

def verify_streamflow_outputs():
    log("=== Verifying HyMAP2 Streamflow Outputs ===")
    results = {}
    for exp in EXPERIMENTS:
        out_dir = os.path.join(P1_DIR, "experiments/NorthMor/routing/output/P3_2019", exp, "ROUTING")
        files = glob.glob(f"{out_dir}/**/*.nc", recursive=True)
        results[exp] = {
            "dir": out_dir,
            "count": len(files),
            "files": sorted(files),
            "status": "OK" if len(files) >= 365 else "INCOMPLETE"
        }
        log(f"Experiment {exp}: {len(files)} streamflow files in {out_dir} (Status: {results[exp]['status']})")
    return results

def main():
    log("Starting 2019 Automated Remapping and Routing Orchestrator...")
    
    # Wait until all 3 DA chains complete
    while True:
        completed = {exp: check_chain_completed(exp) for exp in EXPERIMENTS}
        all_done = all(completed.values())
        
        active = get_active_jobs()
        active_da = {jid: info for jid, info in active.items() if any(exp in info['name'] for exp in ['DA_nocdf', 'DA_cdftrain', 'DA_v2', 'DA_'])}
        
        log(f"DA Status: {completed} | Active DA SLURM jobs: {len(active_da)}")
        
        # Check if any chain failed (no job running, but not completed)
        for exp, done in completed.items():
            if not done:
                is_running = any(exp in info['name'] for info in active_da.values())
                if not is_running:
                    err_file, err_content = check_for_errors(exp)
                    if err_file:
                        log(f"ERROR detected in {exp} ({err_file})! Job is not running.")
        
        if all_done:
            log("SUCCESS: All 3 DA 2019 chains completed successfully through 2020-01-01!")
            break
            
        time.sleep(120)

    # 1. Run remapping
    run_remapping()
    
    # 2. Submit HyMAP jobs
    submitted_ids = submit_hymap_jobs()
    
    # 3. Monitor HyMAP jobs
    monitor_hymap_jobs(submitted_ids)
    
    # 4. Verify streamflow outputs
    streamflow_results = verify_streamflow_outputs()
    
    # 5. Write final status report
    report_path = os.path.join(WORKSPACE, "reports/hymap2_2019_runs_status.md")
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w") as f:
        f.write("# Rapport d'exécution du routage HyMAP 2019\n\n")
        f.write(f"- Date d'achèvement : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## 1. Statut des chaînes DA 2019\n")
        for exp in EXPERIMENTS:
            f.write(f"- `{exp}` : Complétée avec succès (restart 2020-01-01 présent)\n")
        f.write("\n## 2. Jobs HyMAP exécutés\n")
        for exp, jid in submitted_ids.items():
            f.write(f"- `{exp}` : Job SLURM ID `{jid}`\n")
        f.write("\n## 3. Fichiers de débit HyMAP (Streamflow)\n")
        for exp, res in streamflow_results.items():
            f.write(f"- `{exp}` : {res['count']} fichiers .nc dans `{res['dir']}` (Statut : {res['status']})\n")
            if res['count'] > 0:
                f.write(f"  - Premier fichier : `{os.path.basename(res['files'][0])}`\n")
                f.write(f"  - Dernier fichier : `{os.path.basename(res['files'][-1])}`\n")
        f.write("\n*Rappel : Aucun critère ni diagnostic de la Décision 16 n'a été calculé conformément aux consignes.*\n")
    
    log(f"Status report written to {report_path}. Orchestrator finished.")

if __name__ == "__main__":
    main()
