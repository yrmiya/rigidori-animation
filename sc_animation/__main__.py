if __name__ == '__main__':
    import sys
    try:
        from .cli import main
        main()

        exit()

    except KeyboardInterrupt:
        sys.stdout.write('\nAborted\n')
        sys.exit(0)
