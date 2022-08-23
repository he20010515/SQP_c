/*
 * @Author: HeYuwei
 * @Date: 2022-08-18 18:49:16
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-08-19 09:53:17
 * @FilePath: \SQP_c\src\util\hash_table.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include <stdlib.h>
#include <stdio.h>
#include "elog.h"
#include "hash_table.h"
#include "string.h"

#define INIT_TABLE_SIZE 16

/**
 * @description: 申请一个哈希表
 * @param {function} hash_function 哈希函数
 * @param {function} compare_function 比较函数
 * @return {*}
 */
HashTable *HashTable_alloc(unsigned int (*hash_function)(void *),
                           int (*compare_function)(void *, void *))
{

    HashTable *ht = (HashTable *)malloc(sizeof(HashTable));
    if (!ht)
    {
        log_e("alloc HashTable faild");
        return NULL;
    }

    ht->table = malloc(sizeof(HashNode) * INIT_TABLE_SIZE);
    memset(ht->table, 0, INIT_TABLE_SIZE * sizeof(HashNode));
    ht->_buffer_size = INIT_TABLE_SIZE;
    ht->hash_function = hash_function;
    ht->compare_function = compare_function;
    ht->elem_num = 0;

    if (!ht->table)
    {
        log_e("no enough memory for hashtable");
        free(ht);
        ht = NULL;
        return NULL;
    }
    return ht;
}
/**
 * @description: 初始化哈希节点
 * @param {HashNode} *hashnode
 * @return {*}
 */
void HashNode_init(HashNode *hashnode)
{
    hashnode->status = VOID;
    hashnode->key = NULL;
    hashnode->value = NULL;
    hashnode->next = NULL;
}
/**
 * @description: 向哈希表中添加元素
 * @param {HashTable} *hs_table 哈希表
 * @param {void} *key   键
 * @param {void} *value 值
 * @return {*}
 */
int HashTable_insert(HashTable *hs_table, void *key, void *value)
{
    if (!hs_table)
    {
        log_e("error");
        return -1;
    }

    unsigned int i = hs_table->hash_function(key) % hs_table->_buffer_size;
    HashNode *p = hs_table->table[i];
    HashNode *pri = p;

    while (p)
    {
        if (hs_table->compare_function(key, p->key) == 0)
        {
            if (p->status == VOID)
            {
                p->value = value;
                p->status = USED;
                p->key = key;
                break;
            }
        }
        pri = p;
        p = p->next;
    }

    if (!p)
    {
        HashNode *temp = (HashNode *)malloc(sizeof(HashNode));
        if (!temp)
        {
            log_e("no enough memory\n");
            return -1;
        }
        HashNode_init(temp);
        temp->key = key;
        temp->value = value;
        temp->status = USED;

        if (pri == NULL)
        {
            hs_table->table[i] = temp;
        }
        else
        {
            pri->next = temp;
        }
    }
    return 0;
}

void *HashTable_get(HashTable *hs_table, void *key)
{
    if (!hs_table || !key)
    {
        log_e("error");
        return NULL;
    }
    int i = hs_table->hash_function(key) % hs_table->_buffer_size;
    HashNode *temp = hs_table->table[i];
    while (temp)
    {
        if (hs_table->compare_function(temp->key, key) == 1)
        {
            return temp->value;
        }
        temp = temp->next;
    }
    return NULL;
}

void HashTable_free(HashTable *hs_table)
{
    if (hs_table)
    {
        if (hs_table->table)
        {
            for (int i = 0; i < hs_table->_buffer_size; i++)
            {
                HashNode **head = &(hs_table->table[i]);
                HashNode *q = NULL;
                while (*head != NULL)
                {
                    q = *head;
                    *head = q->next;
                    free(q);
                }
            }
            free(hs_table->table);
        }
        free(hs_table);
    }
    return;
}
