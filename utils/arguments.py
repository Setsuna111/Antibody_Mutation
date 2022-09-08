import argparse


def get_common_args():
    parser = argparse.ArgumentParser()
    # the environment setting
    parser.add_argument('--batch_size', type=int, default=4, help='the batch size for once training')
    parser.add_argument('--lr', type=float, default=1e-4, help='the initial learning rate')
    parser.add_argument('--n_epoch', type=int, default=100, help='the number of epoch for training')
    parser.add_argument('--start_epoch', type=int, default=0, help='the number of trained epoch')
    parser.add_argument('--load_checkpoint', type=bool, default=False, help='whether to load checkpoint')
    parser.add_argument('--is_cuda', type=bool, default=False, help='whether to use cuda')

    args = parser.parse_args()
    return args