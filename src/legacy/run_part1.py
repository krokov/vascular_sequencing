import SEVtras
import multiprocessing as mp

def main():
    SEVtras.sEV_recognizer(
        input_path   = "/Users/paulcahill/Outs",
        sample_file  = "/Users/paulcahill/sEV_Run/sample_file",
        out_path     = "/Users/paulcahill/sEV_Run/outputs",
        species      = "Mus",
        dir_origin   = False
        # , threads=1   # uncomment this line if multiprocessing errors continue
    )
    print("✅ Part I finished. Outputs saved in /Users/paulcahill/sEV_Run/outputs")

if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()

