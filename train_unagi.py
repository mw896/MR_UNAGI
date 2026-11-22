import sys
from pathlib import Path

def train(stage_dir, idrem_dir, seed):
    from UNAGI.UNAGI_tool import UNAGI
    stages, idrem = Path(stage_dir).resolve(), Path(idrem_dir).resolve()
    inputs = [stages / f"{i}.h5ad" for i in range(4)]
    for path in (*inputs, idrem / "idrem.jar", idrem / "example_settings.txt"):
        if not path.is_file():
            raise FileNotFoundError(path)
    model = UNAGI(random_seed=seed)
    model.register_species("Human")
    model.setup_data(stages, stage_key="stage", total_stage=4, label_key="marker_driven_celltype")
    if model.input_dim != 30232:
        raise ValueError("Expected 30,232 input genes")
    model.setup_training(task="MM_UNAGI_smallrun", dist="ziln", device="cpu", GPU=False,
                         latent_dim=64, graph_dim=1024, hidden_dim=256, BATCHSIZE=16384,
                         GCN=True, adversarial=True, beta=1, lr=1e-4, lr_dis=5e-4,
                         max_iter=5, epoch_initial=5, epoch_iter=5)
    model.run_UNAGI(idrem_dir=str(idrem), CPO=True)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("Usage: python train_unagi.py STAGE_DIR IDREM_DIR SEED")
    train(sys.argv[1], sys.argv[2], int(sys.argv[3]))
