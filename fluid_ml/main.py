import torch


# some day this will be the entrypoint to a beautiful sciml project yay

def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)


if __name__ == "__main__":
    main()
