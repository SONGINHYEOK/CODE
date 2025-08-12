from celery import shared_task
import subprocess
import os
import tempfile
from pathlib import Path
from django.conf import settings
from .models import ScreeningJob, ScreeningResult
from apps.molecules.services import MoleculeService

@shared_task(bind=True)
def run_unidock_screening(self, job_id):
    """Run UniDock virtual screening"""
    try:
        job = ScreeningJob.objects.get(id=job_id)
        job.status = 'running'
        job.save()
        
        # Create temp directory for UniDock
        with tempfile.TemporaryDirectory() as temp_dir:
            # Prepare input files
            receptor_path = os.path.join(temp_dir, 'receptor.pdbqt')
            ligand_path = os.path.join(temp_dir, 'ligands.sdf')
            output_path = os.path.join(temp_dir, 'output.sdf')
            
            # Copy receptor file
            with open(receptor_path, 'wb') as f:
                f.write(job.receptor_file.read())
            
            # Copy ligand file
            if job.ligand_file:
                with open(ligand_path, 'wb') as f:
                    f.write(job.ligand_file.read())
            
            # Prepare UniDock command
            unidock_cmd = [
                os.path.join(settings.UNIDOCK_PATH, 'unidock'),
                '--receptor', receptor_path,
                '--gpu_batch', ligand_path,
                '--center_x', str(job.batch.center_x),
                '--center_y', str(job.batch.center_y),
                '--center_z', str(job.batch.center_z),
                '--size_x', str(job.batch.size_x),
                '--size_y', str(job.batch.size_y),
                '--size_z', str(job.batch.size_z),
                '--exhaustiveness', str(job.exhaustiveness),
                '--num_modes', str(job.num_modes),
                '--energy_range', str(job.energy_range),
                '--dir', temp_dir,
                '--out', output_path,
            ]
            
            # Add GPU flag if enabled
            if settings.UNIDOCK_GPU_ENABLED:
                unidock_cmd.append('--gpu')
            
            # Run UniDock
            result = subprocess.run(
                unidock_cmd,
                capture_output=True,
                text=True,
                cwd=temp_dir
            )
            
            if result.returncode != 0:
                raise Exception(f"UniDock failed: {result.stderr}")
            
            # Parse results
            parse_unidock_results(job, output_path)
            
            # Save output file
            with open(output_path, 'rb') as f:
                job.result_file.save('results.sdf', f)
            
            job.status = 'completed'
            job.save()
            
    except Exception as e:
        job.status = 'failed'
        job.error_message = str(e)
        job.save()
        raise

def parse_unidock_results(job, output_file):
    """Parse UniDock output and save to database"""
    # This is a simplified parser - adjust based on actual UniDock output format
    from rdkit import Chem
    
    suppl = Chem.SDMolSupplier(output_file)
    
    for mol in suppl:
        if mol is None:
            continue
        
        # Get SMILES
        smiles = Chem.MolToSmiles(mol)
        
        # Get or create molecule
        molecule = MoleculeService.get_or_create_molecule(smiles)
        if not molecule:
            continue
        
        # Get docking score (adjust property name based on UniDock output)
        score = float(mol.GetProp('docking_score') if mol.HasProp('docking_score') else 0)
        
        # Create result
        ScreeningResult.objects.create(
            job=job,
            molecule=molecule,
            docking_score=score,
            is_hit=(score < -7.0)  # Example threshold
        )