#include "DNA.h"
#include "DNA-functions.h"

int IOReadOneOption(FILE* OptionsFile, char* stuk)
{
  char ch, ch2;
  int line;
  char str[DNA_STRINGLENGTH_SPRINTF];

  line = 1;
  while ((ch = getc(OptionsFile)) != EOF)
  {
    if (ch != '\n' && ch != '#')
    {
      while (ch == ' ' || ch == '\t' || ch == '\n')
      {
        ch = getc(OptionsFile);
      }
      /* We now have a non space character */
      memcpy(stuk, &ch, sizeof(char));
      ch2 = getc(OptionsFile);
      ungetc(ch2, OptionsFile);
      if (ch2 != ' ' && ch2 != '\n')
      {
        IOLineGetSkip(str, OptionsFile);
        memcpy(stuk + 1, str, sizeof(str));
      }
      else
      {
        str[0] = '\0';
        stuk[1] = '\0';
      }
      return (line);
    }
    if (ch != '\n')
    {
      IOLineGet(str, OptionsFile);
    }
    line++;
  }

  return (EOF);
}

int IOLineGet(char* ssring, FILE* fp)
{
  char ch;
  int index = 0;

  while ((ch = getc(fp)) != EOF)
  {
    ssring[index++] = ch;
    if (index == DNA_STRINGLENGTH_SPRINTF)
    {
      IOErrorOnScreen(1, "Keyword in file too long!");
    }

    if (ch == '\n')
    {
      ssring[index] = '\0';
      return (index - 1);
    }
  }

  return 0;
}

int IOLineGetSemi(char* ssring, FILE* fp)
{
  char ch;
  int index = 0;

  while ((ch = getc(fp)) != EOF)
  {
    if (ch == ':')
    {
      ssring[index] = '\0';
      return (1);
    }
    else if (ch == '\n')
    {
      ssring[index] = '\0';
      return (0);
    }
    else
    {
      ssring[index++] = ch;
      if (index == DNA_STRINGLENGTH_SPRINTF)
      {
        IOErrorOnScreen(1, "Keyword in file too long!");
      }
    }
  }

  return 0;
}

int IOLineGetSkip(char* ssring, FILE* fp)
{
  char ch;
  int index = 0;

  while ((ch = getc(fp)) != EOF)
  {
    if (ch != ' ' && ch != '\n')
    {
      ssring[index++] = ch;
      if (index == DNA_STRINGLENGTH_SPRINTF)
      {
        IOErrorOnScreen(1, "Keyword in file too long!");
      }
    }
    if ((ch == '\n' || ch == ' ') && index > 0)
    {
      ssring[index] = '\0';
      return (index);
    }
  }

  return 0;
}
