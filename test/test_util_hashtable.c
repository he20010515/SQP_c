/*
 * @Author: HeYuwei
 * @Date: 2022-08-19 09:10:39
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-08-19 09:57:28
 * @FilePath: \SQP_c\test\test_util_hashtable.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include "stdio.h"
#include "hash_table.h"
#include "string.h"
typedef struct commodity
{
    char name[32];
    float price;
} Com;

long long hash(void *key)
{
    int key_len = 32;
    unsigned int hash = 1315423911;
    int i = 0;
    char *str_key = (char *)key;

    for (i = 0; i < key_len; str_key++, i++)
    {
        hash ^= ((hash << 5) + (*str_key) + (hash >> 2));
    }

    return hash;
}
int comparefunction(void *key1, void *key2)
{
    char *_key1, *_key2;
    _key1 = (char *)key1;
    _key2 = (char *)key2;

    return strcmp(key1, key2) == 0;
}

int main(int argc, char const *argv[])
{
    Com coms[50] = {0};
    char buffer[32];
    for (int i = 0; i < 50; i++)
    {
        sprintf(buffer, "coms%d", i);
        strcpy(coms[i].name, buffer);
        coms[i].price = (float)i * i;
    }

    for (int j = 0; j < 3; j++)
    {
        HashTable *table = HashTable_alloc(hash, comparefunction);
        for (int i = 0; i < 50; i++)
        {
            HashTable_insert(table, &(coms[i].name), &(coms[i].price));
        }

        for (int i = 0; i < 50; i++)
        {
            float temp;
            temp = *(float *)HashTable_get(table, coms[i].name);
            printf("key = %s,value=%f,\n", coms[i].name, temp);
        }
        HashTable_free(table);
    }

    return 0;
}