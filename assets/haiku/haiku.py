import random
import re

import numpy as np
from bs4 import BeautifulSoup
import pyscript as ps
from js import window
# import json
from time import perf_counter

# Set up a CORS proxy URL
cors_proxy_url = "https://corsproxy.io/?"
# cors_proxy_url = "http://api.allorigins.win/get?url="

async def get_random_chemrxiv_paper():

    # Send request to ChemRXiv Organic Chemisty website
    url = "https://chemrxiv.org/engage/chemrxiv/category-dashboard/605c72ef153207001f6470d1"

    # Send request to ChemRXiv website via CORS proxy
    print('--> Sending ChemRXiv page request')
    start = perf_counter()

    response = await window.fetch(f"{cors_proxy_url}{url}")
    string = await response.text()

    print(f'--> main page loaded ({perf_counter()-start:.1f} s)')
    print(string[0:100])

    # content = json.loads(string)["contents"]
    content = string
    # with open('main_page.html', 'w') as f:
    #     f.write(html)

    # Parse HTML content using BeautifulSoup
    soup = BeautifulSoup(content, 'html.parser')
    # soup = BeautifulSoup(response.content, 'html.parser')

    # Find all article links on the page
    article_links = soup.find_all('a', class_='ArticleSummary')

    if not article_links:
        return "title", "authors", "abstract"

    # Select a random article link
    random_article_link = random.choice(article_links)

    # Send request to the random article page
    article_url = random_article_link['href']

    print(f'--> Sending ChemRXiv paper page request ({cors_proxy_url}https://chemrxiv.org/{article_url})')
    time = perf_counter()

    # article_response = requests.get(f"{cors_proxy_url}https://chemrxiv.org/"+article_url)
    article_response = await window.fetch(f"{cors_proxy_url}https://chemrxiv.org/"+article_url)
    article_string = await article_response.text()
    # article_content = json.loads(article_string)["contents"]
    article_content = article_string

    print(f'--> article page loaded ({perf_counter()-start:.1f} s)')

    # with open('article_page.html', 'w') as f:
        # f.write(article_response.text)

    # Parse HTML content of the article page
    article_soup = BeautifulSoup(article_content, 'html.parser')

    # Extract title
    title = article_soup.title.text.split('|')[0].rstrip()
    print("--> title is:", title)

    authors = article_soup.select('#main-content > div > div > div.row.LayoutGutters.mx-0.align-center.justify-center.mb-5 > ' +
                                'div > div > div.col-md-8.col-12 > div.article-header > div:nth-child(1) > div:nth-child(2) > ' +
                                'div > div.mt-2.col > ul')[0].text.rstrip()
    # correct for wrong positioning of comma
    authors = authors.replace(' ,', ', ')

    # replace last comma with 'and', unless single author
    authors = ", ".join(authors.split(", ")[:-1]) + ' and ' + authors.split(", ")[-1] if "," in authors else authors

    # get abstract in a text format
    abstract = article_soup.find('div', class_='abstract').text.strip()

    # return the extracted information
    return title, authors, abstract

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
            return (random.choice((
                "Sometimes even I,\nlike the very best, struggle.\nLet me try again.",
                "It's never easy\nto find all you seek in life.\nI should try again.",
                "Give me another\nchance at composing haikus.\nYou won't regret it.",
                "Roses are not red,\nnor violetes are really blue.\nDid not find haikus.",
                "Such essential words.\nIt's less easy with less text.\nLet me try again.",
                )), 2)

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
    input_div = ps.document.querySelector("#text")
    string = input_div.value.replace("\n", " ")

    contiguous = ps.document.querySelector("#checkbox").checked

    output_div = ps.document.querySelector("#output")
    output_div.innerText, exit_code = get_haiku(string, contiguous=contiguous)

    output_div_extra = ps.document.querySelector("#output_extra")
    output_div_extra.innerText = ("",
                                    "[no haiku present]",
                                    "[max iterations reached, retry]")[exit_code]

async def scrape_chemrxiv(event):

    input_div = ps.document.querySelector("#text")
    title, authors, input_div.innerText = await get_random_chemrxiv_paper()

    string = input_div.value.replace("\n", " ")

    checkbox = ps.document.querySelector("#checkbox")
    checkbox.checked = False

    output_div = ps.document.querySelector("#output")
    output_div.innerText, exit_code = get_haiku(string, contiguous=False, maxiter=1E4)

    output_div_extra = ps.document.querySelector("#output_extra")
    output_div_extra.innerText = ("",
                                    "[no haiku present]",
                                    "[max iterations reached, retry]")[exit_code]
    
    credits_div = ps.document.querySelector("#paper_credits")
    credits_div.innerText = f"from \"{title}\" by {authors}"