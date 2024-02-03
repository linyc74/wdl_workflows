import time
import argparse
import subprocess


REQUIRED = [
    {
        'keys': ['--host'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'host name of the container registry (e.g. *.azurecr.io)',
        }
    },
    {
        'keys': ['--repo'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'container repository name',
        }
    },
    {
        'keys': ['--version'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'container version',
        }
    },
    {
        'keys': ['--ubuntu-version'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'ubuntu version',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        time_stamp = time.strftime('%Y-%m-%d-%H-%M')
        cmd = f'docker build -t {args.host}/{args.repo}:{args.version}_{args.ubuntu_version}_{time_stamp} .'
        print(f'\n{cmd}\n', flush=True)
        subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    EntryPoint().main()
