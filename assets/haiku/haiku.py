import numpy as np
import re
from pyscript import document

# masterpiece single liner to count syllables from 200_success
# https://codereview.stackexchange.com/questions/224177/python-function-to-count-syllables-in-a-word
def count_syllables(word):
    return len(
        re.findall('(?!e$)[aeiouy]+', word, re.I) +
        re.findall('^[^aeiouy]*e$', word, re.I)
    )

def capitalize(verse):
    pieces = verse.split()
    pieces[0] = pieces[0].capitalize()
    return " ".join(pieces)

bad_endings = (
    'and',
    'of',
    'the',
)

def get_haiku(string, contiguous=True, maxiter=5E3):
    '''
    return: haiku, exit code (0: normal, 1: no haiku present, 2: maxiter reached)

    '''

    # initialize haiku container
    haiku = ["", "", ""]

    # keep track of the position of first words in verses
    # so that we ensure no repetitions
    first_words_indices = set()

    # split line after opening or closing round or square brackets and after spaces
    splitted = [word for word in re.split("\[|\]|\(|\)| ", string) if word not in ("", ".", ",")]

    # create a randomly shuffled index iterator
    index_iterator = iter(np.random.permutation(np.arange(len(splitted))))

    # set the initial index cursor to None
    index = None

    # iteration counter to prevent being stuck should
    # anything go wrong (because of course it will)
    iteration = 0

    while True:
        
        # check if we have the correct number of syllables
        if [sum(count_syllables(w) for w in v.split()) for v in haiku] == [5, 7, 5]:

            # make sure that the last word is not a weird one - if not, we got it
            if haiku[2].split()[-1].rstrip('.').rstrip(',') not in bad_endings:
                break

            # if it is a bad ending, delete everything (if contiguous)...
            if contiguous:
                haiku = ["", "", ""]

            # ...or just the last verse (if not) and try to find a better one
            else:
                haiku[2] = ""
            
        # safety break out of main while loop
        iteration += 1
        if iteration == maxiter:
            return ("Sometimes even I,\nlike the very best, struggle.\nLet me try again.", 2)

        try:

            # get the verse index
            verse_index = next(v for v, verse in enumerate(haiku) if verse == "")

            # get a new random index for a word
            # if we are starting from scratch or 
            # for every new verse in non-contiguous mode
            if not contiguous or haiku[0] == "" or index == None:
                while index is None:
                    
                    # get a new index from the iterator
                    temp_index = next(index_iterator)

                    # make sure we have not started there already
                    if temp_index not in first_words_indices:

                        # make sure this word is not a number
                        if re.search("\d", splitted[temp_index]) is None:
                            index = temp_index
                            first_words_indices.add(temp_index)

            # make a line out of that word (check for numbers again)
            while re.search("\d", splitted[index]) is not None:
                index += 1
            line = splitted[index]

            # add words until we have enough (or too many) syllables
            while sum(count_syllables(w) for w in line.split()) < [5, 7, 5][verse_index]:

                # move cursor by one position
                index += 1

                # make sure to avoid numbers
                while re.search("\d", splitted[index]) is not None:
                    index += 1

                # add word to verse
                line += " " + splitted[index]
            
            # if the verse has the right length, save it and move the cursor
            if sum(count_syllables(w) for w in line.split()) == [5, 7, 5][verse_index]:
                haiku[verse_index] = line

                # if we retain the verse move the cursor
                #  so that we do not repeat words
                index += 1

            # if we failed to get a verse of the right length and we
            # want a continuous one, reset the haiku and start fresh
            elif contiguous:
                haiku = ["", "", ""]

            # if non-contiguous, reset iterator, since we can use
            # any starting point again
            else:
                index_iterator = iter(np.random.permutation(np.arange(len(splitted))))

        # if we run out of words at the end of the string at any point, restart
        # from verse (non-contiguous) or entire haiku (contiuguous) and go on.
        # We stop trying only if we started from every possible point.
        except IndexError:
            index = None

        except StopIteration:
            break

        # Give up trying if we did what we could
        if len(first_words_indices) == len(splitted):
            break

    if [sum(count_syllables(w) for w in line.split()) for line in haiku] == [5, 7, 5]:

        # capitalize the first word of the first verse
        haiku[0] = capitalize(haiku[0])

        # remove trailing comma, if present
        haiku[2] = haiku[2][:-1] if haiku[2][-1] == "," else haiku[2]

        # Add a full stop at the end if not there already
        haiku[2] += "." if haiku[2][-1] != "." else ""

        return ("\n".join(haiku), 0)
    
    return ("No haiku was found.\nPlease, try adding some more text.\nAmusement awaits.", 1)

def main(event):
    input_text = document.querySelector("#text")
    string = input_text.value.replace("\n", " ")

    contiguous = document.querySelector("#checkbox").checked

    output_div = document.querySelector("#output")
    output_div.innerText, exit_code = get_haiku(string, contiguous=contiguous)

    output_div_extra = document.querySelector("#output_extra")
    output_div_extra.innerText = ("",
                                    "[no haiku present]",
                                    "[max iterations reached, retry]")[exit_code]