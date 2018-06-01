import sys

def progress(total, progress, epoch, status):
    """
    Displays or updates a console progress bar.

    Original source: https://stackoverflow.com/a/15860757/1391441
    """
    barLength = 20
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(barLength * progress))
    text = "\r[Epoch {epoch}]: [{bar}] {percentage}% {status}".format(
        epoch = epoch,
        bar = "#" * block + "-" * (barLength - block), 
        percentage = round(progress * 100, 0),
        status = status)
    sys.stdout.write(text)
    if progress < 1:
        sys.stdout.flush()