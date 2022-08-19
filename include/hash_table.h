/*
 * @Author: HeYuwei
 * @Date: 2022-08-19 09:09:04
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-08-19 09:40:40
 * @FilePath: \SQP_c\include\hash_table.h
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
typedef enum HashNodeStatus
{
    VOID = 0,
    USED
} HashNodeStatus;

typedef struct hash_node
{
    void *key;
    void *value;
    HashNodeStatus status;
    struct hash_node *next;

} HashNode;

typedef struct hash_table
{
    long long _buffer_size;
    long long elem_num;
    long long (*hash_function)(void *);
    int (*compare_function)(void *, void *); //比较函数,若相等返回1
    HashNode **table;
} HashTable;

HashTable *HashTable_alloc(long long (*hash_function)(void *),
                           int (*compare_function)(void *, void *));

int HashTable_insert(HashTable *hs_table, void *key, void *value);
void *HashTable_get(HashTable *hs_table, void *key);
void HashTable_free(HashTable *hs_table);